use std::f64::consts::PI;

use faer::linalg::triangular_solve::solve_upper_triangular_in_place;
use faer::mat::{AsMatMut, AsMatRef};
use faer::{mat, perm, unzip, zip, Col, ColRef, Mat, MatRef, Par, Side};
use indicatif::ParallelProgressIterator;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::constants::{J, MU_B, SCALAR_J};
use crate::utils::*;
use crate::{Coupling, MagneticField, C64};

/// Minimum energy that isn't just set for zero (in meV)
const ZERO_ENERGY_TOL: f64 = 1e-12;
const SINGULAR_FTOL: f64 = 1e-7;
const SINGULAR_CTOL: C64 = C64::new(1e-7, 0.);
const J2PI: C64 = C64::new(0., 2. * PI);

/// The result of a single-Q spinwave calculation.
///
/// Fields:
/// - `energies`: The energies of the spinwave modes at this q-vector.
/// - `sab`: Optional spin-spin correlation tensors S'^{alpha, beta} for each mode.
/// - `intensities`: The neutron scattering intensities (S_perp) for each mode.
pub struct SpinwaveResult {
    pub energies: Vec<C64>,
    pub sab: Option<Vec<Mat<C64>>>,
    pub intensities: Vec<f64>,
    pub wavefunctions: Option<Mat<C64>>
}

/// The q-independent components of the calculation.
/// Fields:
/// - `C`: The q-independent component of the Hamiltonian matrix.
/// - `z`: The z-components (R_x + i R_y) of the rotation matrices.
/// - `sab_blocks`: The q-independent parts of Y, Z, V, W, blocks from computing Sab
/// - `spin_coefficients`: The spin coefficients matrix.
/// - `Az`: Optional Zeeman term for the `A` matrix if a magnetic field is provided.
struct QIndependentComponents {
    C: Mat<C64>,
    z: Vec<Col<C64>>,
    sab_blocks: Mat<Mat<C64>>,
    spin_coefficients: Mat<C64>,
    Az: Option<Vec<C64>>,
}

/// Components required for the rotating frame calculations (Toth & Lake eq 39)
/// Fields:
/// - `km`: The propagation vector
/// - `nx`: The cross-product matrix of the perpendicular vector in Rodrigues' formula
/// - `R1`: The rotation matrix for the +/-km components
/// - `R2`: The rotation matrix for the k=0 component
struct RotatingFrameComponents {
    km: Col<f64>,
    nx: Mat<C64>,
    R1: Mat<C64>,
    R2: Mat<C64>,
}

/// Return values for hermitian/non-hermitian branches - eigenvalues and transformation matrix T
/// Fields:
/// - `eigvals`: The eigenvalues of the Hamiltonian
/// - `T`: The transformation matrix to compute the spin correlations (derived from eigenvectors)
struct HamiltonianResults {
    eigvals: Col<C64>,
    T: Mat<C64>,
}

/// Calculate the q-independent components of the calculation.
///
/// # Parameters
/// - `rotations`: Rotation matrices for each site.
/// - `magnitudes`: Magnitudes of the moments at each site.
/// - `couplings`: Slice of couplings between sites.
/// - `field`: Optional external magnetic field.
///
/// # Returns
/// The q-independent components of the calculation.
#[inline(always)]
fn calc_q_independent(
    rotations: Vec<MatRef<C64>>,
    magnitudes: Vec<f64>,
    couplings: &[&Coupling],
    field: Option<MagneticField>,
) -> QIndependentComponents {
    let n_sites = rotations.len();

    // decompose rotation matrices
    // in the notation of Petit (2011)
    // etas[i] is the direction of the i'th moment in Cartesian coordinates
    let (z, etas) = get_rotation_components(rotations);

    // make spin coefficients array
    // so spin_coefficients[i, j] = sqrt(S_i S_j) / 2
    let root_mags = Col::<C64>::from_iter(magnitudes.iter().map(|x| C64::from((0.5 * x).sqrt())));
    let spin_coefficients = root_mags.clone() * root_mags.transpose();

    // create matrix C of Hamiltonian which is q-independent
    let mut C = Mat::<C64>::zeros(n_sites, n_sites);
    for c in couplings {
        C[(c.index2, c.index2)] += spin_coefficients[(c.index2, c.index2)]
            * (etas[c.index1].transpose() * c.matrix.as_ref() * etas[c.index2]);
    }
    C *= 2.;

    // Adds a small delta to diagonal to ensure we don't have an exact degeneracies
    C.diagonal_mut()
        .column_vector_mut()
        .iter_mut()
        .enumerate()
        .for_each(|(i, x)| *x -= C64::from((i as f64) * ZERO_ENERGY_TOL));

    // if an external magnetic field is provided, calculate the Az matrix
    // which is added to the diagonal of A in the Hamiltonian
    let Az: Option<Vec<C64>> = match field {
        Some(f) => Some(
            f.g_tensors
                .into_iter()
                .enumerate()
                .map(|(i, t)| -0.5 * MU_B * f.vector.transpose() * t * etas[i])
                .collect(),
        ),
        None => None,
    };

    // We store S' as a 3x3 matrix of matrices indexed by alpha, beta
    // where each element is an 2 * n_sites x 2 * n_sites matrix
    // later on we will sum over the diagonal elements to get S'^alpha,beta(k, omega)
    let mut sab_blocks = Mat::<Mat<C64>>::full(3, 3, Mat::<C64>::new());

    // note one can show:
    // Y*[alpha, beta] = Y[beta, alpha]
    // Z*[alpha, beta] = V[beta, alpha]
    // V*[alpha, beta] = Z[beta, alpha]
    // W*[alpha, beta] = W[beta, alpha]
    // thus we only need to calculate one triangle and can fill in the rest by conjugation
    for alpha in 0..3 {
        let z_alphas = Col::<C64>::from_iter(z.iter().map(|zi| zi[alpha]));
        for beta in 0..=alpha {
            let z_betas = Col::<C64>::from_iter(z.iter().map(|zi| zi[beta]));

            // construct the four blocks V, W, Y, Z;
            // note that before we include the phase factors,
            // V is conj(Z) and W is conj(Y)
            let Yab = z_alphas.clone() * z_betas.adjoint();
            let Zab = z_alphas.clone() * z_betas.transpose();
            let Vab = Zab.conjugate().to_owned();
            let Wab = Yab.conjugate().to_owned();

            sab_blocks[(alpha, beta)] = block_matrix(&Yab, &Zab, &Vab, &Wab);

            if beta < alpha {
                let Yba = Yab.adjoint().to_owned();
                let Zba = Vab.adjoint().to_owned();
                let Vba = Zab.adjoint().to_owned();
                let Wba = Wab.adjoint().to_owned();

                sab_blocks[(beta, alpha)] = block_matrix(&Yba, &Zba, &Vba, &Wba);
            }
        }
    }

    QIndependentComponents {
        C,
        z,
        sab_blocks,
        spin_coefficients,
        Az,
    }
}

/// Calculate the spin Hamiltonian matrix at a particular q-vector
///
/// # Parameters
/// - `q`: The q-vector for which to calculate the Hamiltonian.
/// - `q_independent_components`: The q-independent component of the Hamiltonian.
/// - `n_sites`: The number of sites in the system.
/// - `couplings`: Slice of couplings between sites.
///
/// # Returns
/// The spin Hamiltonian matrix h(q)
#[inline(always)]
fn calc_spin_hamiltonian(
    q: Col<f64>,
    q_independent_components: &QIndependentComponents,
    n_sites: usize,
    couplings: &[&Coupling],
) -> Mat<C64> {
    // unpack q-independent components
    let C = &q_independent_components.C;
    let z = &q_independent_components.z;
    let spin_coefficients = &q_independent_components.spin_coefficients;
    let Az = &q_independent_components.Az;

    // create A and B matrices for the Hamiltonian

    let mut A = Mat::<C64>::zeros(n_sites, n_sites);
    let mut B = Mat::<C64>::zeros(n_sites, n_sites);

    for c in couplings {
        let phase_factor = ((2. * J * PI) * (q.transpose() * c.inter_site_vector.as_ref())).exp();
        let (i, j) = (c.index1, c.index2);

        // contributions to A and B from this coupling
        A[(i, j)] += (z[i].transpose() * c.matrix.as_ref() * z[j].conjugate()) * phase_factor;
        B[(i, j)] += (z[i].transpose() * c.matrix.as_ref() * z[j].as_ref()) * phase_factor;
    }

    A = component_mul(&A, spin_coefficients);
    B = component_mul(&B, spin_coefficients);

    // slightly convoluted way to add to the diagonal of A because adding a Diag to a Mat
    // isn't implemented by `faer` yet (missed out, apparently!)
    if let Some(Az_vals) = Az {
        for i in 0..n_sites {
            A[(i, i)] += Az_vals[i];
        }
    }

    let A_minus_C: Mat<C64> = A.clone() - C;
    let A_conj_minus_C: Mat<C64> = A.adjoint() - C;
    let B_adj = B.adjoint().to_owned();

    block_matrix(&A_minus_C, &B, &B_adj, &A_conj_minus_C)
}

/// Calculate energies and correlation S'^{alpha, beta} for a set of q-vectors.
///
/// # Parameters
/// - `rotations`: Rotation matrices for each site.
/// - `magnitudes`: Magnitudes of the moments at each site.
/// - `q_vectors`: A vector of q-vectors for which to calculate energies and intensities.
/// - `couplings`: Slice of couplings between sites.
/// - `positions`: The relative positions of each site within the unit cell.
/// - `rlu_to_cart`: Optional matrix to transform q to Cartesian for intensity calculation.
/// - `field`: Optional external magnetic field.
/// - `rotating_frame`: Optional pair of vectors (kvec, nperp) defining a rotating frame
/// - `save_Sab`: Boolean indicating if full Sab matrices are returned (default: False)
///
///# Returns
/// A tuple containing:
/// - A vector of vectors, where each inner vector contains the energies for the corresponding
///   q-vector.
/// - A vector of vectors, where each inner vector contains S'^{alpha, beta} matrices for each
///   eigenvalue at the corresponding q-vector.
pub fn calc_spinwave(
    rotations: Vec<MatRef<C64>>,
    magnitudes: Vec<f64>,
    q_vectors: Vec<Vec<f64>>,
    in_couplings: Vec<&Coupling>,
    positions: Vec<ColRef<f64>>,
    rlu_to_cart: Option<MatRef<f64>>,
    field: Option<MagneticField>,
    rotating_frame: Option<Vec<ColRef<f64>>>,
    save_Sab: bool,
    save_wavefunctions: bool,
) -> Vec<SpinwaveResult> {
    let n_sites = rotations.len();
    let n_q = q_vectors.len() as u64;

    let mut couplings = in_couplings;
    // Need a second vector with mutable elements from which we can borrow due to Pyo3 limits.
    let new_couplings: Vec<Coupling>;
    let rotating_components: Option<RotatingFrameComponents> = match rotating_frame {
        Some(rot_comps) => {
            let km: Col<f64> = rot_comps[0].to_owned();
            let n = Col::<C64>::from_iter(rot_comps[1].iter().map(|f| C64::new(*f, 0.)));
            // Computes the rotation matrices in Toth & Lake eq 39
            let nx = mat![
                [C64::ZERO, -n[2], n[1]],
                [n[2], C64::ZERO, -n[0]],
                [-n[1], n[0], C64::ZERO]
            ];
            let R2 = n.clone() * n.transpose();
            let R1 = (Mat::<C64>::identity(3, 3) - (nx.as_ref() * SCALAR_J) - R2.as_ref()) / 2.;
            new_couplings = Vec::from_iter(couplings.iter().map(|c| {
                let phi = 2. * PI * km.transpose() * c.inter_site_vector.as_ref();
                // R is the Rodrigues rotation matrix
                let R = Mat::<C64>::identity(3, 3) * phi.cos()
                    + &nx * phi.sin()
                    + (1. - phi.cos()) * &R2;
                Coupling {
                    index1: c.index1,
                    index2: c.index2,
                    matrix: (c.matrix.as_ref() * R.as_ref() + R.as_ref() * c.matrix.as_ref()) / 2.,
                    inter_site_vector: c.inter_site_vector.clone(),
                }
            }));
            couplings = Vec::from_iter(new_couplings.iter());
            Some(RotatingFrameComponents { km, nx, R1, R2 })
        }
        _ => None,
    };

    let QIndependentComponents = calc_q_independent(rotations, magnitudes, &couplings, field);

    if rotating_components.is_some() {
        // Caculate energies and intensities for a triplet of q-vectors in a rotating frame
        q_vectors
            .into_par_iter()
            .progress_count(n_q)
            .map(|q| {
                let q_vec = Col::from_iter(q);
                (-1..=1)
                    .into_par_iter()
                    .map(|tri_id| {
                        spinwave_single_q(
                            q_vec.clone(),
                            &QIndependentComponents,
                            n_sites,
                            &couplings,
                            &positions,
                            rlu_to_cart,
                            &rotating_components,
                            save_Sab,
                            save_wavefunctions,
                            tri_id as f64,
                        )
                    })
                    .reduce_with(
                        |mut triplet_results: SpinwaveResult, result: SpinwaveResult| {
                            triplet_results.energies.extend(result.energies);
                            if save_Sab {
                                triplet_results
                                    .sab
                                    .as_mut()
                                    .unwrap()
                                    .extend(result.sab.unwrap());
                            }
                            triplet_results.intensities.extend(result.intensities);
                            triplet_results
                        },
                    )
                    .unwrap()
            })
            .collect()
    } else {
        q_vectors
            .into_par_iter()
            .progress_count(n_q)
            .map(|q| {
                spinwave_single_q(
                    Col::from_iter(q),
                    &QIndependentComponents,
                    n_sites,
                    &couplings,
                    &positions,
                    rlu_to_cart,
                    &rotating_components,
                    save_Sab,
                    save_wavefunctions,
                    0.,
                )
            })
            .collect()
    }
}

/// Tries to calculate the sqrt_hamiltonian using cholesky-type decompositions
///
/// # Parameters
/// - `hamiltonian`: The Hamiltonian matrix
///
/// Routine will return error if hamiltonian is not positive semi-definite
fn calc_sqrt_hamiltonian(hamiltonian: MatRef<C64>) -> Result<Mat<C64>, String> {
    // take square root of Hamiltonian using Cholesky if possible; if this fails,
    // use the LDL (Bunch-Kaufmann) decomposition instead and take sqrt(H) = L * sqrt(D)
    if let Ok(chol) = hamiltonian.clone().llt(Side::Lower) {
        Ok(chol.L().to_owned())
    } else if let Ok(ldl) = hamiltonian.ldlt(Side::Lower) {
        let D = ldl.D().column_vector();
        if D.iter().any(|v| v.re < -SINGULAR_FTOL) {
            return Err(String::from("LDLT: Hamiltonian not positive-semi-definite"));
        }
        let mut sqrt_d = Col::<C64>::zeros(hamiltonian.nrows());
        zip!(&mut sqrt_d, D).for_each(|unzip!(sqd, v)| *sqd = v.sqrt());
        Ok(ldl.L() * sqrt_d.as_diagonal())
    } else {
        /* TODO - work out difference between Python and Rust w.r.t. twin_example.py
        let ldl = hamiltonian.lblt(Side::Lower);
        let l = ldl.L();
        let d = ldl.B_diag().column_vector();
        let s = ldl.B_subdiag().column_vector();
        // The Bunch-Kaufmann algorithm works for indefinite matrices.
        // But if the hamiltonian constructed above is not positive-semi-definite
        // then g*ham will have imaginary eigenvalues and we need to use a different
        // algorithm. We can check for semi-positive-definiteness here by checking if
        // any off-diagonal element of the block-diagonal matrix are not zero
        // or if any diagonal element is negative. I've not found any proof that the
        // diagonal blocks are equivalent to eigenvalues but this post:
        // https://mathoverflow.net/questions/84420 suggest it could be.
        if s.iter().any(|v| v.re.abs() > SINGULAR_FTOL) || d.iter().any(|v| v.re < -SINGULAR_FTOL) { */
            return Err(String::from("LBLT: Hamiltonian not positive-semi-definite"));
       /* } 
        let mut sqrt_d = Col::<C64>::zeros(hamiltonian.nrows());
        zip!(&mut sqrt_d, d).for_each(|unzip!(sqd, v)| *sqd = v.sqrt());
        // need to apply permutations: in Python scipy does this for you
        let mut shm = Mat::<C64>::zeros(l.nrows(), l.ncols());
        perm::permute_cols(
            shm.as_mat_mut(),
            (ldl.P().inverse() * l * sqrt_d.as_diagonal()).as_mat_ref(),
            ldl.P().inverse(),
        );
        Ok(shm) */
    }
}

/// Calculate the transition eigenvalues and matrix T using the hermitian algorithm
///
/// # Parameters
/// - `hamiltonian`: The Hamiltonian matrix
/// - `n_sites`: The number of sites in the system.
fn solve_ham_hermitian(mut sqrt_hamiltonian: Mat<C64>, n_sites: usize) -> HamiltonianResults {
    let mut shc: Mat<C64> = sqrt_hamiltonian.clone();
    let mut negative_half = shc.submatrix_mut(n_sites, 0, n_sites, 2 * n_sites);
    negative_half *= -1.;

    let eigendecomp = (sqrt_hamiltonian.adjoint() * shc)
        .self_adjoint_eigen(Side::Lower)
        .expect("Could not calculate eigendecomposition of the Hamiltonian.");

    let eigvals: Col<C64> = eigendecomp.S().column_vector().to_owned();

    // we reverse the columns of U to match the nonincreasing order of eigenvalues in sqrt_E
    let eigvecs: MatRef<C64> = eigendecomp.U().reverse_cols();

    // calculate transformation matrix for spin-spin correlation function
    // this is T = K^-1 U sqrt(E) where E is the diagonal 2 * n_sites matrix of eigenvalues
    // where the first n_sites entries are sqrt(eigval) and the remaining are sqrt(-eigval)
    // for the eigenvalues of the Hamiltonian
    // these should be in nonincreasing order but it doesn't matter as the array comes out
    // the same once we apply the commutation sign flips
    let mut sqrt_E = eigvals.clone();
    sqrt_E.iter_mut().for_each(|x| {
        *x = match *x {
            x if x.re < -ZERO_ENERGY_TOL => (-x).sqrt(),
            x if x.re < ZERO_ENERGY_TOL => C64::ZERO,
            _ => x.sqrt(),
        }
    });

    // Check if sqrt_hamiltonian is singular and if so add delta to diagonal to avoid this
    // As sqrt_hamiltonian is triangular, it is singular if any diagonal element is zero
    if sqrt_hamiltonian
        .diagonal()
        .column_vector()
        .iter()
        .any(|v| v.re.abs() < SINGULAR_FTOL)
    {
        sqrt_hamiltonian
            .diagonal_mut()
            .column_vector_mut()
            .iter_mut()
            .for_each(|x| *x += SINGULAR_CTOL);
    }

    // instead of inverting K and calculating T = K^-1 U sqrt(E),
    // it's faster and more stable to solve the linear system K T = U sqrt(E)
    // note the `faer` solver is in-place so calculates it directly on the variable `T`
    // (the input T is initially the righthand side of the equation U sqrt(E))
    let mut T = eigvecs * sqrt_E.as_diagonal();
    solve_upper_triangular_in_place(sqrt_hamiltonian.adjoint().as_ref(), T.as_mut(), Par::Seq);
    HamiltonianResults { eigvals, T }
}

/// Calculate the eigenvalues and transition matrix T using the non-hermitian algorithm
///
/// # Parameters
/// - `hamiltonian`: The Hamiltonian matrix
/// - `n_sites`: The number of sites in the system.
fn solve_ham_nonherm(hamiltonian: Mat<C64>, n_sites: usize) -> HamiltonianResults {
    let mut hc: Mat<C64> = hamiltonian.clone();
    let mut negative_half = hc.submatrix_mut(n_sites, 0, n_sites, 2 * n_sites);
    negative_half *= -1.;
    let eigendecomp = hc
        .eigen()
        .expect("Could not calculate eigendecomposition of the Hamiltonian.");
    let eigvals: Col<C64> = eigendecomp.S().column_vector().to_owned();
    let eigvecs: MatRef<C64> = eigendecomp.U();
    let mut T = eigvecs.cloned();
    // Computes gComm * V * gComm
    let mut botleft = T.submatrix_mut(n_sites, 0, n_sites, n_sites);
    botleft *= -1.;
    let mut topright = T.submatrix_mut(0, n_sites, n_sites, n_sites);
    topright *= -1.;
    let mut M = (T.adjoint() * eigvecs)
        .diagonal()
        .column_vector()
        .to_owned();
    M.iter_mut()
        .for_each(|m| *m = (C64::ONE / (*m + SINGULAR_CTOL)).sqrt());
    T = eigvecs * M.as_diagonal();
    HamiltonianResults { eigvals, T }
}

/// Calculate energies and intensities for a single q-vector.
///
/// # Parameters
/// - `q`: The q-vector for which to calculate energies and intensities.
/// - `q_independent_components`: The q-independent components for the calculation.
/// - `n_sites`: The number of sites in the system.
/// - `couplings`: Slice of couplings between sites.
/// - `positions`: The relative positions of each site within the unit cell.
/// - `rlu_to_cart`: Optional matrix to transform q to Cartesian for intensity calculation.
/// - `rotating_components`: Option with fields needed for rotating frame calculation.
/// - `save_sab`: Whether to save the 3x3 Sab tensors or not (default: False).
/// - `tri_id`: For optional rotating frame calculation, whether this is -k, 0 or +k
///
///# Returns
/// A tuple containing:
/// - A vector containing the energies for the given q-vector.
/// - A vector of S'^{alpha, beta} matrices for each eigenvalue at the given q-vector.
fn spinwave_single_q(
    mut q: Col<f64>,
    q_independent_components: &QIndependentComponents,
    n_sites: usize,
    couplings: &[&Coupling],
    positions: &[ColRef<f64>],
    rlu_to_cart: Option<MatRef<f64>>,
    rotating_components: &Option<RotatingFrameComponents>,
    save_Sab: bool,
    save_wavefunctions: bool,
    tri_id: f64,
) -> SpinwaveResult {
    let sab_blocks = &q_independent_components.sab_blocks;
    let spin_coefficients = &q_independent_components.spin_coefficients;

    if let Some(rotcomp) = rotating_components {
        q += tri_id * rotcomp.km.as_ref();
    }

    let hamiltonian =
        calc_spin_hamiltonian(q.clone(), q_independent_components, n_sites, couplings);
    let ham_results = if let Ok(sqrt_hamiltonian) = calc_sqrt_hamiltonian(hamiltonian.as_ref()) {
        solve_ham_hermitian(sqrt_hamiltonian, n_sites)
    } else {
        // If the hamiltonian is not positive semi-definite we use the non-hermitian algorithm
        solve_ham_nonherm(hamiltonian, n_sites)
    };
    let eigvals = &ham_results.eigvals;
    let T = &ham_results.T;

    // calculate block matrices [ Y Z ; V W ] for S'^alpha, beta
    //
    // The calculation includes the phase factors exp(i q (r_i - r_j)) for each
    // pair of sites i, j; to calculate these efficiently we calculate
    // exp(i q r_i) for each site i and then the outer product of this with its conjugate
    // gives us the full matrix of phase factors.
    let phase_factors = Col::<C64>::from_iter(
        positions
            .iter()
            .map(|r_i| (J2PI * (q.transpose() * r_i)).exp()),
    );
    let phase_factors_matrix = phase_factors.clone() * phase_factors.adjoint();

    let coeffblk = component_mul(&(2. * spin_coefficients), &phase_factors_matrix);
    let coefficients = block_matrix(&coeffblk, &coeffblk, &coeffblk, &coeffblk);

    // Apply transformation matrix to S'^alpha,beta block matrices T*[VW;YZ]T
    // and then we just take the diagonal elements as that's all we need for
    // S'^alpha,beta(k, omega) at each eigenvalue
    // We do the division required by 2 * n_sites here to save doing it later
    let block_diags = Mat::<Col<C64>>::from_fn(3, 3, |alpha, beta| -> Col<C64> {
        let mat =
            T.adjoint() * component_mul(&sab_blocks[(alpha, beta)], &coefficients) * T.as_ref();
        mat.diagonal().column_vector().to_owned() / (2 * n_sites) as f64
    });

    // now create S' for each eigenvalue (the only places where there are non-zero intensities)
    let mut Sab: Vec<Mat<C64>> = (0..2 * n_sites)
        .map(|i| {
            // each element of S' over alpha, beta is created from an index over 2 * n_sites
            Mat::<C64>::from_fn(3, 3, |alpha, beta| -> C64 {
                block_diags[(alpha, beta)].as_ref()[i]
            })
        })
        .collect();

    // For rotating frame calculation, apply rotation transformations to Sab
    if let Some(rotcomp) = rotating_components {
        // Convert back into lab frame (eq 37)
        let R2I = &rotcomp.R2 - Mat::<C64>::identity(3, 3);
        let R22I = (2. * &rotcomp.R2) - Mat::<C64>::identity(3, 3);
        Sab.iter_mut().for_each(|m| {
            *m = 0.5
                * (&*m - (&rotcomp.nx * &*m * &rotcomp.nx)
                    + (&R2I * &*m * &rotcomp.R2)
                    + (&rotcomp.R2 * &*m * &R22I))
        });
        // Apply the rotation transformation (eq 40)
        match tri_id {
            n if n < 0. => Sab
                .iter_mut()
                .for_each(|m| *m = &*m * rotcomp.R1.conjugate()),
            n if n > 0. => Sab.iter_mut().for_each(|m| *m = &*m * &rotcomp.R1),
            0. => Sab.iter_mut().for_each(|m| *m = &*m * &rotcomp.R2),
            _ => panic!(),
        }
        // Convert q back for qperp calculation
        q -= tri_id * rotcomp.km.as_ref();
    }

    // gets the conversion from r.l.u. to Cartesian, for Sperp we need Q in Cartesians
    let qcart: Col<f64> = match rlu_to_cart {
        Some(m) => (q.transpose() * m).transpose().to_owned(),
        None => q,
    };

    // and finally, calculate the perpendicular component of Sab
    let intensities = {
        let mut norm_q: Col<C64> = (qcart.as_ref() / qcart.norm_l2())
            .iter()
            .map(C64::from)
            .collect();
        if norm_q.has_nan() {
            norm_q = Col::<C64>::from_iter([C64::ZERO, C64::ZERO, C64::ZERO]);
        }
        let perp_factor = Mat::<C64>::identity(3, 3) - (norm_q.as_ref() * norm_q.adjoint());
        Sab.iter()
            .map(|Sab_w| -> f64 { component_mul(Sab_w, &perp_factor).sum().re })
            .collect()
    };

    let sab_out = match save_Sab {
            true => Some(Sab),
            false => None
    };

    let wavefunctions_out = match save_wavefunctions {
        true => Some(T.clone()),
        false => None
    };

    SpinwaveResult {
            energies: eigvals.iter().copied().collect(),
            sab: sab_out,
            intensities: intensities,
            wavefunctions: wavefunctions_out
    }
}
