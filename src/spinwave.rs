use std::f64::consts::PI;

use faer::{linalg::solvers::Solve, unzip, zip, Col, ColRef, Mat, MatRef, Side};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::constants::{J, MU_B};
use crate::utils::{component_mul, get_rotation_components, make_block_hamiltonian};
use crate::{Coupling, MagneticField, C64};

/// Result of a spinwave calculation.
pub struct SpinwaveResult {
    pub energies: Vec<f64>,
    pub correlation: Option<Vec<Mat<C64>>>,
}

/// Calculate the q-independent components of the calculation.
fn calc_q_independent(
    rotations: Vec<MatRef<C64>>,
    magnitudes: Vec<f64>,
    couplings: &Vec<&Coupling>,
    field: Option<MagneticField>,
) -> (Mat<C64>, Vec<Col<C64>>, Mat<C64>, Option<Vec<C64>>) {
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
    (C, z, spin_coefficients, Az)
}

/// Calculate the square root of the Hamiltonian.
#[inline(always)]
fn calc_sqrt_hamiltonian(
    q: Col<f64>,
    C: &Mat<C64>,
    n_sites: usize,
    z: &[Col<C64>],
    spin_coefficients: &Mat<C64>,
    couplings: &[&Coupling],
    Az: &Option<Vec<C64>>,
) -> Mat<C64> {
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

    let hamiltonian: Mat<C64> = make_block_hamiltonian(A, B, C);

    // take square root of Hamiltonian using Cholesky if possible; if this fails,
    // use the LDL (Bunch-Kaufmann) decomposition instead and take sqrt(H) = L * sqrt(D)
    let sqrt_hamiltonian = {
        if let Ok(chol) = hamiltonian.clone().llt(Side::Lower) {
            chol.L().to_owned()
        } else {
            let ldl = hamiltonian.lblt(Side::Lower);
            let l = ldl.L();
            let d = ldl.B_diag().column_vector(); // we're ignoring off-diagonals... this may be
                                                  // dangerous

            // we use the zip and unzip to map over d and allocate to sqrt_d
            let mut sqrt_d = Col::<C64>::zeros(d.nrows());
            zip!(&mut sqrt_d, d).for_each(|unzip!(sqd, v)| *sqd = v.sqrt());

            // need to apply permutations: in Python scipy does this for you
            ldl.P().inverse() * l * sqrt_d.as_diagonal()
        }
    };
    sqrt_hamiltonian
}

pub fn calc_energies(
    rotations: Vec<MatRef<C64>>,
    magnitudes: Vec<f64>,
    q_vectors: Vec<Vec<f64>>,
    couplings: Vec<&Coupling>,
    field: Option<MagneticField>,
) -> Vec<Vec<f64>> {
    let n_sites = rotations.len();

    let (C, z, spin_coefficients, Az) =
        calc_q_independent(rotations, magnitudes, &couplings, field);

    // now perform the calculation for each q-vector in parallel
    q_vectors
        .into_par_iter()
        .map(|q| {
            energies_single_q(
                Col::from_iter(q),
                &C,
                n_sites,
                &z,
                &spin_coefficients,
                &couplings,
                &Az,
            )
        })
        .collect()
}

/// Calculate energies (eigenvalues of the Hamiltonian) for a single q-value
fn energies_single_q(
    q: Col<f64>,
    C: &Mat<C64>,
    n_sites: usize,
    z: &[Col<C64>],
    spin_coefficients: &Mat<C64>,
    couplings: &[&Coupling],
    Az: &Option<Vec<C64>>,
) -> Vec<f64> {
    let sqrt_hamiltonian =
        calc_sqrt_hamiltonian(q, C, n_sites, z, spin_coefficients, couplings, Az);
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

/// Calculate the block matrices for S'^alpha, beta
/// That is, the matrix [ Y Z ; V W ] for each alpha, beta pair
/// note that the phase factor is NOT applied here; this is done later
fn calc_sab_blocks(
    n_sites: usize,
    z: &[Col<C64>],
    q: Col<f64>,
    spin_coefficients: &Mat<C64>,
    positions: &[ColRef<f64>],
) -> Mat<Mat<C64>> {
    let phase_factors = Col::<C64>::from_iter(
        positions
            .iter()
            .map(|r_i| (J * (q.transpose() * r_i)).exp()),
    );
    let phase_factors_matrix = phase_factors.clone() * phase_factors.conjugate().transpose();

    let coefficients = component_mul(&spin_coefficients, &phase_factors_matrix);

    // We store S' as a 3x3 matrix of matrices indexed by alpha, beta
    // where each element is an 2 * n_sites x 2 * n_sites matrix
    // later on we will sum over the diagonal elements to get S'^alpha,beta(k, omega)
    Mat::<Mat<C64>>::from_fn(3, 3, |alpha, beta| -> Mat<C64> {
        let mut matrix = Mat::<C64>::zeros(2 * n_sites, 2 * n_sites);
        let z_alphas = Col::<C64>::from_iter(z.iter().map(|zi| zi[alpha]));
        let z_betas = Col::<C64>::from_iter(z.iter().map(|zi| zi[beta]));

        // construct the four blocks V, W, Y, Z;
        // note that while we're just dealing in the z components,
        // V is conj(Z) and W is conj(Y)
        let Y = z_alphas.clone() * z_betas.conjugate().transpose();
        let Z = z_alphas * z_betas.transpose();
        let V = Z.conjugate().to_owned();
        let W = Y.conjugate().to_owned();

        component_mul(&Y, &coefficients);
        component_mul(&Z, &coefficients);
        component_mul(&V, &coefficients);
        component_mul(&W, &coefficients);

        // copy the blocks into the matrix
        matrix.submatrix_mut(0, 0, n_sites, n_sites).copy_from(&Y);
        matrix
            .submatrix_mut(0, n_sites, n_sites, n_sites)
            .copy_from(&Z);
        matrix
            .submatrix_mut(n_sites, 0, n_sites, n_sites)
            .copy_from(&V);
        matrix
            .submatrix_mut(n_sites, n_sites, n_sites, n_sites)
            .copy_from(&W);

        matrix
    })
}

pub fn calc_spinwave(
    rotations: Vec<MatRef<C64>>,
    magnitudes: Vec<f64>,
    q_vectors: Vec<Vec<f64>>,
    couplings: Vec<&Coupling>,
    positions: Vec<ColRef<f64>>,
    field: Option<MagneticField>,
) -> Vec<SpinwaveResult> {
    let n_sites = rotations.len();

    let (C, z, spin_coefficients, Az) =
        calc_q_independent(rotations, magnitudes, &couplings, field);

    q_vectors
        .into_par_iter()
        .map(|q| {
            spinwave_single_q(
                Col::from_iter(q),
                &C,
                n_sites,
                &z,
                &spin_coefficients,
                &couplings,
                &positions,
                &Az,
            )
        })
        .collect()
}

/// Calculate energies and intensities for a single q-point.
fn spinwave_single_q(
    q: Col<f64>,
    C: &Mat<C64>,
    n_sites: usize,
    z: &[Col<C64>],
    spin_coefficients: &Mat<C64>,
    couplings: &[&Coupling],
    positions: &[ColRef<f64>],
    Az: &Option<Vec<C64>>,
) -> SpinwaveResult {
    let sqrt_hamiltonian =
        calc_sqrt_hamiltonian(q.clone(), C, n_sites, z, spin_coefficients, couplings, Az);
    let mut shc: Mat<C64> = sqrt_hamiltonian.clone();
    let mut negative_half = shc.submatrix_mut(n_sites, 0, n_sites, 2 * n_sites);
    negative_half *= -1.;

    let sab_blocks = calc_sab_blocks(n_sites, z, q, spin_coefficients, positions);

    let eigendecomp = (sqrt_hamiltonian.adjoint() * shc)
        .self_adjoint_eigen(Side::Lower)
        .expect("Could not calculate eigendecomposition of the Hamiltonian.");

    let eigvals: ColRef<C64> = eigendecomp.S().column_vector();

    // we reverse the rows of U to match the nonincreasing order of eigenvalues in sqrt_E
    let eigvecs: MatRef<C64> = eigendecomp.U().reverse_rows();

    // calculate transformation matrix for spin-spin correlation function
    // this is T = K^-1 U sqrt(E) where E is the diagonal 2 * n_sites matrix of eigenvalues
    // where the first n_sites entries are sqrt(eigval) and the remaining are sqrt(-eigval)
    // for the eigenvalues of the Hamiltonian
    // we use `reverse_rows()` to sort eigvals in nonincreasing order (default is nondecreasing order)
    // as Toth & Lake (2015) assumes eigenvalues are nonincreasing
    let mut sqrt_E = eigvals.reverse_rows().to_owned();
    let mut negative_half = sqrt_E.subrows_mut(n_sites, n_sites);
    negative_half *= -1.;
    sqrt_E.iter_mut().for_each(|x| *x = x.sqrt());

    // instead of inverting K and calculating T = K^-1 U sqrt(E),
    // it's faster and more stable to solve the linear system K T = U sqrt(E)
    // note the `faer` solver is in-place so calculates it directly on the variable `T`
    // (the input T is initially the righthand side of the equation U sqrt(E))
    let mut T = eigvecs * sqrt_E.as_diagonal();
    sqrt_hamiltonian.partial_piv_lu().solve_in_place(T.as_mut());

    // Apply transformation matrix to S'^alpha,beta block matrices T*[VW;YZ]T
    // and then we just take the diagonal elements as that's all we need for
    // S'^alpha,beta(k, omega) at each eigenvalue
    let block_diags = Mat::<Col<C64>>::from_fn(3, 3, |alpha, beta| -> Col<C64> {
        let mat = T.adjoint() * sab_blocks[(alpha, beta)].as_ref() * T.as_ref();
        mat.diagonal().column_vector().to_owned()
    });

    // now create S' for each eigenvalue (the only places where there are non-zero intensities)
    let Sab: Vec<Mat<C64>> = (0..n_sites)
        .map(|i| {
            // each element of S' over alpha, beta is created from an index over 2 * n_sites
            Mat::<C64>::from_fn(3, 3, |alpha, beta| -> C64 {
                let diag: ColRef<C64> = block_diags[(alpha, beta)].as_ref();
                diag[i] / (2 * n_sites) as f64
            })
        })
        .chain((0..n_sites).map(|i| {
            // negative energy modes
            Mat::<C64>::from_fn(3, 3, |alpha, beta| -> C64 {
                let diag: ColRef<C64> = block_diags[(alpha, beta)].as_ref();
                diag[n_sites - i] / (2 * n_sites) as f64
            })
        }))
        .collect();

    SpinwaveResult {
        energies: eigvals.iter().map(|x| x.re).collect(),
        correlation: Some(Sab),
    }
}
