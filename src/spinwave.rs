use std::f64::consts::PI;

use faer::linalg::triangular_solve::solve_upper_triangular_in_place;
use faer::{unzip, zip, Col, ColRef, Mat, MatRef, Par, Side, perm};
use faer::mat::{AsMatRef, AsMatMut};
use indicatif::ParallelProgressIterator;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::constants::{J, MU_B};
use crate::utils::*;
use crate::{Coupling, MagneticField, C64};

/// Minimum energy that isn't just set for zero (in meV)
const ZERO_ENERGY_TOL: f64 = 1e-12;

/// The result of a single-Q spinwave calculation.
///
/// Fields:
/// - `energies`: The energies of the spinwave modes at this q-vector.
/// - `sab`: Optional spin-spin correlation tensors S'^{alpha, beta} for each mode.
/// - `intensities`: The neutron scattering intensities (S_perp) for each mode.
pub struct SpinwaveResult {
    pub energies: Vec<f64>,
    pub sab: Option<Vec<Mat<C64>>>,
    pub intensities: Vec<f64>,
}

/// The q-independent components of the calculation.
/// Fields:
/// - `C`: The q-independent component of the Hamiltonian matrix.
/// - `z`: The z-components (R_x + i R_y) of the rotation matrices.
/// - `spin_coefficients`: The spin coefficients matrix.
/// - `Az`: Optional Zeeman term for the `A` matrix if a magnetic field is provided.
struct QIndependentComponents {
    C: Mat<C64>,
    z: Vec<Col<C64>>,
    spin_coefficients: Mat<C64>,
    Az: Option<Vec<C64>>,
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

    QIndependentComponents {
        C,
        z,
        spin_coefficients,
        Az,
    }
}

/// Calculate the square root of the Hamiltonian.
///
/// # Parameters
/// - `q`: The q-vector for which to calculate the Hamiltonian.
/// - `C`: The q-independent component of the Hamiltonian matrix.
/// - `n_sites`: The number of sites in the system.
/// - `z`: The z-components (R_x + i R_y) of the rotation matrices.
/// - `spin_coefficients`: The spin coefficients matrix.
/// - `couplings`: Slice of couplings between sites.
/// - `Az`: Optional Zeeman term for the `A` matrix if a magnetic field is provided.
///
/// # Returns
/// The "square root" of the Hamiltonian matrix (K such that K K* = h(q)).
#[inline(always)]
fn calc_sqrt_hamiltonian(
    q: Col<f64>,
    q_indepepdent_components: &QIndependentComponents,
    n_sites: usize,
    couplings: &[&Coupling],
) -> Mat<C64> {
    // unpack q-independent components
    let C = &q_indepepdent_components.C;
    let z = &q_indepepdent_components.z;
    let spin_coefficients = &q_indepepdent_components.spin_coefficients;
    let Az = &q_indepepdent_components.Az;

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

    // Adds a small delta to diagonal to ensure we don't have an exact degeneracies
    let diag_delta = Mat::<C64>::from_fn(A.nrows() * 2, A.ncols() * 2,
        |i,j| if i == j { C64::from((i as f64)*1e-12) } else { C64::ZERO } );
    let hamiltonian: Mat<C64> = block_matrix(&A_minus_C, &B, &B_adj, &A_conj_minus_C) + diag_delta;

    // take square root of Hamiltonian using Cholesky if possible; if this fails,
    // use the LDL (Bunch-Kaufmann) decomposition instead and take sqrt(H) = L * sqrt(D)
    if let Ok(chol) = hamiltonian.clone().llt(Side::Lower) {
        chol.L().to_owned()
    } else {
        if let Ok(ldl) = hamiltonian.ldlt(Side::Lower) {
            let mut sqrt_d = Col::<C64>::zeros(hamiltonian.nrows());
            zip!(&mut sqrt_d, ldl.D().column_vector()).for_each(|unzip!(sqd, v)| *sqd = v.sqrt());
            ldl.L() * sqrt_d.as_diagonal()
        } else {
            let ldl = hamiltonian.lblt(Side::Lower);
            let l = ldl.L();
            let d = ldl.B_diag().column_vector(); // we're ignoring off-diagonals... this may be
                                                  // dangerous
            let mut sqrt_d = Col::<C64>::zeros(hamiltonian.nrows());
            zip!(&mut sqrt_d, d).for_each(|unzip!(sqd, v)| *sqd = v.sqrt());
            // need to apply permutations: in Python scipy does this for you
            let mut shm = Mat::<C64>::zeros(l.nrows(), l.ncols());
            perm::permute_cols(shm.as_mat_mut(), (ldl.P().inverse() * l * sqrt_d.as_diagonal()).as_mat_ref(), ldl.P().inverse());
            shm
        }} 
}

/// Calculate energies (eigenvalues of the Hamiltonian) for a set of q-vectors.
///
/// # Parameters
/// - `rotations`: Rotation matrices for each site.
/// - `magnitudes`: Magnitudes of the moments at each site.
/// - `q_vectors`: A vector of q-vectors for which to calculate energies.
/// - `couplings`: Slice of couplings between sites.
/// - `field`: Optional external magnetic field.
///
/// # Returns
/// A vector of vectors, where each inner vector contains the energies for the corresponding q-vector.
pub fn calc_energies(
    rotations: Vec<MatRef<C64>>,
    magnitudes: Vec<f64>,
    q_vectors: Vec<Vec<f64>>,
    couplings: Vec<&Coupling>,
    field: Option<MagneticField>,
) -> Vec<Vec<f64>> {
    let n_sites = rotations.len();
    let n_q = q_vectors.len() as u64;

    let q_independent_components = calc_q_independent(rotations, magnitudes, &couplings, field);

    // now perform the calculation for each q-vector in parallel
    q_vectors
        .into_par_iter()
        .progress_count(n_q)
        .map(|q| {
            energies_single_q(
                Col::from_iter(q),
                &q_independent_components,
                n_sites,
                &couplings,
            )
        })
        .collect()
}

/// Calculate energies (eigenvalues of the Hamiltonian) for a single q-value.
///
/// # Parameters
/// - `q`: The q-vector for which to calculate energies.
/// - `C`: The q-independent component of the Hamiltonian matrix.
/// - `n_sites`: The number of sites in the system.
/// - `z`: The z-components (R_x + i R_y) of the rotation matrices.
/// - `spin_coefficients`: The spin coefficients matrix.
/// - `couplings`: Slice of couplings between sites.
/// - `Az`: Optional Zeeman term for the `A` matrix if a magnetic field is provided.
///
///# Returns
/// A vector containing the energies for the given q-vector.
fn energies_single_q(
    q: Col<f64>,
    q_independent_components: &QIndependentComponents,
    n_sites: usize,
    couplings: &[&Coupling],
) -> Vec<f64> {
    let sqrt_hamiltonian = calc_sqrt_hamiltonian(q, q_independent_components, n_sites, couplings);
    // 'shc' is "square root of Hamiltonian with commutation"`
    // We need to enforce the bosonic commutation properties, we do this
    // by finding the 'square root' of the matrix (i.e. finding K such that KK^dagger = H)
    // and then negating the second half.
    //
    // In matrix form we do
    //
    //     M = K^dagger g K
    //
    // where g is a diagonal matrix of length 2n, with the first n entries being 1, and the
    // remaining entries being -1.
    // We do this by just multiplying the >n_sites rows of shc to get g*K
    let mut shc: Mat<C64> = sqrt_hamiltonian.clone();
    let mut negative_half = shc.submatrix_mut(n_sites, 0, n_sites, 2 * n_sites);
    negative_half *= -1.;

    (sqrt_hamiltonian.adjoint() * shc)
        .self_adjoint_eigenvalues(Side::Lower)
        .expect("Could not calculate eigendecomposition of the Hamiltonian.")
}

/// Calculate energies and correlation S'^{alpha, beta} for a set of q-vectors.
///
/// # Parameters
/// - `rotations`: Rotation matrices for each site.
/// - `magnitudes`: Magnitudes of the moments at each site.
/// - `q_vectors`: A vector of q-vectors for which to calculate energies and intensities.
/// - `couplings`: Slice of couplings between sites.
/// - `positions`: The relative positions of each site within the unit cell.
/// - `field`: Optional external magnetic field.
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
    couplings: Vec<&Coupling>,
    positions: Vec<ColRef<f64>>,
    rlu_to_cart: MatRef<f64>,
    field: Option<MagneticField>,
    save_Sab: bool,
) -> Vec<SpinwaveResult> {
    let n_sites = rotations.len();
    let n_q = q_vectors.len() as u64;

    let QIndependentComponents = calc_q_independent(rotations, magnitudes, &couplings, field);

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
                save_Sab,
            )
        })
        .collect()
}

/// Calculate energies and intensities for a single q-vector.
///
/// # Parameters
/// - `q`: The q-vector for which to calculate energies and intensities.
/// - `q_independent_components`: The q-independent components for the calculation.
/// - `n_sites`: The number of sites in the system.
/// - `couplings`: Slice of couplings between sites.
/// - `positions`: The relative positions of each site within the unit cell.
///
///# Returns
/// A tuple containing:
/// - A vector containing the energies for the given q-vector.
/// - A vector of S'^{alpha, beta} matrices for each eigenvalue at the given q-vector.
fn spinwave_single_q(
    q: Col<f64>,
    q_independent_components: &QIndependentComponents,
    n_sites: usize,
    couplings: &[&Coupling],
    positions: &[ColRef<f64>],
    rlu_to_cart: MatRef<f64>,
    save_Sab: bool,
) -> SpinwaveResult {
    let z = &q_independent_components.z;
    let spin_coefficients = &q_independent_components.spin_coefficients;

    let mut sqrt_hamiltonian =
        calc_sqrt_hamiltonian(q.clone(), q_independent_components, n_sites, couplings);
    let mut shc: Mat<C64> = sqrt_hamiltonian.clone();
    let mut negative_half = shc.submatrix_mut(n_sites, 0, n_sites, 2 * n_sites);
    negative_half *= -1.;

    // calculate block matrices [ Y Z ; V W ] for S'^alpha, beta
    //
    // The calculation includes the phase factors exp(i q (r_i - r_j)) for each
    // pair of sites i, j; to calculate these efficiently we calculate
    // exp(i q r_i) for each site i and then the outer product of this with its conjugate
    // gives us the full matrix of phase factors.
    let J2PI = 2. * J * PI;
    let phase_factors = Col::<C64>::from_iter(
        positions
            .iter()
            .map(|r_i| (J2PI * (q.transpose() * r_i)).exp()),
    );
    let phase_factors_matrix = phase_factors.clone() * phase_factors.adjoint();

    let coefficients = component_mul(&(2. * spin_coefficients), &phase_factors_matrix);

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
            let mut Yab = z_alphas.clone() * z_betas.adjoint();
            let mut Zab = z_alphas.clone() * z_betas.transpose();
            let mut Vab = Zab.conjugate().to_owned();
            let mut Wab = Yab.conjugate().to_owned();

            Yab = component_mul(&Yab, &coefficients);
            Zab = component_mul(&Zab, &coefficients);
            Vab = component_mul(&Vab, &coefficients);
            Wab = component_mul(&Wab, &coefficients);

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

    let eigendecomp = (sqrt_hamiltonian.adjoint() * shc).self_adjoint_eigen(Side::Lower)
        .expect("Could not calculate eigendecomposition of the Hamiltonian.");

    let eigvals: ColRef<C64> = eigendecomp.S().column_vector();

    // we reverse the columns of U to match the nonincreasing order of eigenvalues in sqrt_E
    let eigvecs: MatRef<C64> = eigendecomp.U().reverse_cols();

    // calculate transformation matrix for spin-spin correlation function
    // this is T = K^-1 U sqrt(E) where E is the diagonal 2 * n_sites matrix of eigenvalues
    // where the first n_sites entries are sqrt(eigval) and the remaining are sqrt(-eigval)
    // for the eigenvalues of the Hamiltonian
    // these should be in nonincreasing order but it doesn't matter as the array comes out
    // the same once we apply the commutation sign flips
    let mut sqrt_E = eigvals.to_owned();
    let mut negative_half = sqrt_E.subrows_mut(0, n_sites);
    negative_half *= -1.;
    sqrt_E.iter_mut().for_each(|x| {
        *x = match *x {
            x if x.re < ZERO_ENERGY_TOL => C64::ZERO,
            _ => x.sqrt(),
        }
    });

    // instead of inverting K and calculating T = K^-1 U sqrt(E),
    // it's faster and more stable to solve the linear system K T = U sqrt(E)
    // note the `faer` solver is in-place so calculates it directly on the variable `T`
    // (the input T is initially the righthand side of the equation U sqrt(E))
    let mut T = eigvecs * sqrt_E.as_diagonal();
    solve_upper_triangular_in_place(sqrt_hamiltonian.adjoint().as_ref(), T.as_mut(), Par::Seq);

    // T is NaN if sqrt_hamiltonian is singular; add a delta to the diagonal to avoid this
    if T.has_nan() {
        sqrt_hamiltonian.diagonal_mut().column_vector_mut().iter_mut().for_each(|x| *x += C64::from(1e-7) );
        T = eigvecs * sqrt_E.as_diagonal();
        solve_upper_triangular_in_place(sqrt_hamiltonian.adjoint().as_ref(), T.as_mut(), Par::Seq);
    }

    // Apply transformation matrix to S'^alpha,beta block matrices T*[VW;YZ]T
    // and then we just take the diagonal elements as that's all we need for
    // S'^alpha,beta(k, omega) at each eigenvalue
    // We do the division required by 2 * n_sites here to save doing it later
    let block_diags = Mat::<Col<C64>>::from_fn(3, 3, |alpha, beta| -> Col<C64> {
        let mat = T.adjoint() * sab_blocks[(alpha, beta)].as_ref() * T.as_ref();
        mat.diagonal().column_vector().to_owned() / (2 * n_sites) as f64
    });

    // now create S' for each eigenvalue (the only places where there are non-zero intensities)
    let Sab: Vec<Mat<C64>> = (0..2 * n_sites)
        .map(|i| {
            // each element of S' over alpha, beta is created from an index over 2 * n_sites
            Mat::<C64>::from_fn(3, 3, |alpha, beta| -> C64 {
                block_diags[(alpha, beta)].as_ref()[i]
            })
        })
        .collect();

    // gets the conversion from r.l.u. to Cartesian, for Sperp we need Q in Cartesians
    let qcart: Col<f64> = (q.transpose() * rlu_to_cart).transpose().to_owned();
    
    // and finally, calculate the perpendicular component of Sab
    let intensities = {
        let mut norm_q: Col<C64> = (qcart.as_ref() / qcart.norm_l2()).iter().map(C64::from).collect();
        if norm_q.has_nan() {
            norm_q = Col::<C64>::from_iter([C64::ZERO, C64::ZERO, C64::ZERO]);
        }
        let perp_factor = Mat::<C64>::identity(3, 3) - (norm_q.as_ref() * norm_q.adjoint());
        Sab.iter()
            .map(|Sab_w| -> f64 { component_mul(Sab_w, &perp_factor).sum().re })
            .collect()
    };

    match save_Sab {
        true => SpinwaveResult {
            energies: eigvals.iter().map(|x| x.re).collect(),
            sab: Some(Sab),
            intensities,
        },
        false => SpinwaveResult {
            energies: eigvals.iter().map(|x| x.re).collect(),
            sab: None,
            intensities,
        },
    }
}
