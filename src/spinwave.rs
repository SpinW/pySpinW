use std::f64::consts::PI;

use faer::{linalg::solvers::DenseSolveCore, unzip, zip, Col, ColRef, Mat, MatRef, Side};
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::constants::{J, MU_B};
use crate::utils::{component_mul, get_rotation_components, make_block_hamiltonian};
use crate::{Coupling, MagneticField, C64};

pub struct SpinwaveResult {
    pub energies: Vec<f64>,
    pub correlation: Option<Vec<Mat<C64>>>,
}

/// Generate a spinwave calculation function for a given single-q calculation function.
///
/// This macro contains the q-independent setup for the spinwave calculation.
macro_rules! gen_spinwave_function {
    ($name:ident, $single_q_fn:ident) => {
        #[allow(non_snake_case)]
        pub fn $name(
            rotations: Vec<MatRef<C64>>,
            magnitudes: Vec<f64>,
            q_vectors: Vec<Vec<f64>>,
            couplings: Vec<&Coupling>,
            field: Option<MagneticField>,
        ) -> Vec<SpinwaveResult> {
            let n_sites = rotations.len();

            // decompose rotation matrices
            // in the notation of Petit (2011)
            // etas[i] is the direction of the i'th moment in Cartesian coordinates
            let (z, etas) = get_rotation_components(rotations);

            // make spin coefficients array
            // so spin_coefficients[i, j] = sqrt(S_i S_j) / 2
            let root_mags =
                Col::<C64>::from_iter(magnitudes.iter().map(|x| C64::from((0.5 * x).sqrt())));
            // we use this `let binding` to create a longer lived value to avoid some compiler complaints
            let binding = root_mags.clone() * root_mags.transpose();
            let spin_coefficients = binding.transpose();

            // create matrix C of Hamiltonian which is q-independent
            let mut C = Mat::<C64>::zeros(n_sites, n_sites);
            for c in &couplings {
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

            // now perform the calculation for each q-vector in parallel
            q_vectors
                .into_par_iter()
                .map(|q| {
                    $single_q_fn(
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
    };
}

gen_spinwave_function!(calc_energies, energies_single_q);
gen_spinwave_function!(calc_spinwave, spinwave_single_q);

/// Calculate the square root of the Hamiltonian,
/// and optionally the block matrix [YZ;VW] for S^alpha,beta
fn calc_sqrt_hamiltonian(
    q: Col<f64>,
    C: &Mat<C64>,
    n_sites: usize,
    z: &[Col<C64>],
    spin_coefficients: &MatRef<C64>,
    couplings: &[&Coupling],
    Az: &Option<Vec<C64>>,
    calc_Sab: bool,
) -> (Mat<C64>, Option<Mat<Mat<C64>>>) {
    // create A and B matrices for the Hamiltonian
    // as well as the workspace for the dynamical correlation function
    // S'^alpha,beta(k, omega)

    let mut A = Mat::<C64>::zeros(n_sites, n_sites);
    let mut B = Mat::<C64>::zeros(n_sites, n_sites);

    // We store S' as a 3x3 matrix of matrices indexed by alpha, beta
    // where each element is an 2 * n_sites x 2 * n_sites matrix
    // later on we will sum over the diagonal elements to get S'^alpha,beta(k, omega)
    let mut S_prime_mats = match calc_Sab {
        true => Some(Mat::<Mat<C64>>::full(
            3,
            3,
            Mat::<C64>::zeros(2 * n_sites, 2 * n_sites),
        )),
        false => None,
    };

    for c in couplings {
        let phase_factor = ((2. * J * PI) * (q.transpose() * c.inter_site_vector.as_ref())).exp();
        let (i, j) = (c.index1, c.index2);

        // contributions to A and B from this coupling
        A[(i, j)] += (z[i].transpose() * c.matrix.as_ref() * z[j].conjugate()) * phase_factor;
        B[(i, j)] += (z[i].transpose() * c.matrix.as_ref() * z[j].as_ref()) * phase_factor;

        if let Some(ref mut spm) = S_prime_mats {
            // contributions to S' from this coupling
            let common_term = spin_coefficients[(i, j)] * phase_factor;
            for alpha in 0..3 {
                for beta in 0..3 {
                    spm[(alpha, beta)][(i, j)] += common_term * (z[i][alpha] * z[j][beta].conj());
                    spm[(alpha, beta)][(i + n_sites, j + n_sites)] +=
                        common_term * (z[i][alpha].conj() * z[j][beta]);
                    spm[(alpha, beta)][(i, j + n_sites)] +=
                        common_term * (z[i][alpha] * z[j][beta]);
                    spm[(alpha, beta)][(i + n_sites, j)] +=
                        common_term * (z[i][alpha].conj() * z[j][beta].conj());
                }
            }
        }
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
    (sqrt_hamiltonian, S_prime_mats)
}

/// Calculate energies (eigenvalues of the Hamiltonian) for a single q-value
fn energies_single_q(
    q: Col<f64>,
    C: &Mat<C64>,
    n_sites: usize,
    z: &[Col<C64>],
    spin_coefficients: &MatRef<C64>,
    couplings: &[&Coupling],
    Az: &Option<Vec<C64>>,
) -> SpinwaveResult {
    let (sqrt_hamiltonian, _) =
        calc_sqrt_hamiltonian(q, C, n_sites, z, spin_coefficients, couplings, Az, false);
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

    SpinwaveResult {
        energies: (sqrt_hamiltonian.adjoint() * shc)
            .self_adjoint_eigenvalues(Side::Lower)
            .expect("Could not calculate eigendecomposition of the Hamiltonian."),
        correlation: None,
    }
}

/// Calculate energies and intensities for a single q-point.
fn spinwave_single_q(
    q: Col<f64>,
    C: &Mat<C64>,
    n_sites: usize,
    z: &[Col<C64>],
    spin_coefficients: &MatRef<C64>,
    couplings: &[&Coupling],
    Az: &Option<Vec<C64>>,
) -> SpinwaveResult {
    let (sqrt_hamiltonian, S_prime_mats) = calc_sqrt_hamiltonian(
        q.clone(),
        C,
        n_sites,
        z,
        spin_coefficients,
        couplings,
        Az,
        true,
    );
    let sab = S_prime_mats.unwrap();

    let mut shc: Mat<C64> = sqrt_hamiltonian.clone();
    let mut negative_half = shc.submatrix_mut(n_sites, 0, n_sites, 2 * n_sites);
    negative_half *= -1.;

    let eigendecomp = (sqrt_hamiltonian.adjoint() * shc)
        .self_adjoint_eigen(Side::Lower)
        .expect("Could not calculate eigendecomposition of the Hamiltonian.");

    let eigvals: ColRef<C64> = eigendecomp.S().column_vector();
    let eigvecs: MatRef<C64> = eigendecomp.U();

    // calculate transformation matrix for spin-spin correlation function
    // this is T = K^-1 * U * sqrt(E) where E is the diagonal 2 * n_sites matrix of eigenvalues
    // where the first n_sites entries are sqrt(eigval) and the remaining are sqrt(-eigval)
    // for the eigenvalues of the Hamiltonian
    let inv_K = sqrt_hamiltonian.partial_piv_lu().inverse();
    let mut sqrt_E = eigvals.to_owned();
    let mut negative_half = sqrt_E.subrows_mut(n_sites, n_sites);
    negative_half *= -1.;
    sqrt_E.iter_mut().for_each(|x| *x = x.sqrt());

    let T: Mat<C64> = inv_K * eigvecs * sqrt_E.as_diagonal();
    // Apply transformation matrix to S'^alpha,beta block matrices T*[VW;YZ]T
    // and then we just take the diagonal elements as that's all we need to calculate
    // S'^alpha,beta(k, omega) at each eigenvalue
    let block_diags = Mat::<Col<C64>>::from_fn(3, 3, |alpha, beta| -> Col<C64> {
        let mat = T.adjoint() * sab[(alpha, beta)].as_ref() * T.as_ref();
        mat.diagonal().column_vector().to_owned()
    });

    // now create S' for each eigenvalue (the only places where there are non-zero intensities)
    let Sab: Vec<Mat<C64>> = (0..n_sites)
        .map(|i| {
            // each element of S' over alpha, beta is created from an index over 2 * n_sites
            Mat::<C64>::from_fn(3, 3, |alpha, beta| -> C64 {
                let diag: ColRef<C64> = block_diags[(alpha, beta)].as_ref();
                diag[i] + 2. * diag[i + n_sites]
            })
        })
        .collect();

    SpinwaveResult {
        energies: eigvals.iter().map(|x| x.re).collect(),
        correlation: Some(Sab),
    }
}
