use std::f64::consts::PI;
use std::iter::zip;

use faer::{concat, unzip, zip, Col, Mat, MatRef, Scale, Side};
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::{Coupling, C64};

// for convenience
// `Scale` is the faer scalar type used for matrix-scalar multiplication
static J: C64 = C64::new(0., 1.);
static SCALAR_J: Scale<C64> = Scale(J);

/// Run the main calculation step for a spinwave calculation.
#[allow(non_snake_case)]
pub fn calc_spinwave(
    rotations: Vec<MatRef<C64>>,
    magnitudes: Vec<f64>,
    q_vectors: Vec<Col<f64>>,
    couplings: Vec<&Coupling>,
) -> Vec<Vec<f64>> {
    let n_sites = rotations.len();

    // decompose rotation matrices
    // in the notation of Petit (2011)
    // eta[i] is the direction of the i'th moment in Cartesian coordinates
    let z: Vec<Col<C64>> = zip(
        _get_rotation_component(&rotations, 0),
        _get_rotation_component(&rotations, 1),
    )
    .map(|(r1, r2)| r1 + (r2 * SCALAR_J))
    .collect();
    let etas = _get_rotation_component(&rotations, 2);

    // make spin coefficients array
    // so spin_coefficients[i, j] = sqrt(S_i S_j) / 2
    let root_mags = Col::<C64>::from_iter(magnitudes.iter().map(|x| C64::from((0.5 * x).sqrt())));
    // we use this `let binding` to create a longer lived value to avoid some compiler complaints
    let binding = root_mags.clone() * root_mags.transpose();
    let spin_coefficients = binding.transpose();

    // create matrix C of Hamiltonian which is q-independent
    let mut C = Mat::<C64>::zeros(n_sites, n_sites);
    for c in &couplings {
        C[(c.index2, c.index2)] += spin_coefficients[(c.index2, c.index2)]
            * (etas[c.index1].transpose() * c.matrix.as_ref() * etas[c.index2].as_ref());
    }
    C *= 2.;

    q_vectors
        .into_par_iter()
        .map(|q| _spinwave_single_q(q, &C, n_sites, &z, &spin_coefficients, &couplings))
        .collect()
}

/// Calculate spectra for a single q-value.
#[allow(non_snake_case)]
fn _spinwave_single_q(
    q: Col<f64>,
    C: &Mat<C64>,
    n_sites: usize,
    z: &[Col<C64>],
    spin_coefficients: &MatRef<C64>,
    couplings: &Vec<&Coupling>,
) -> Vec<f64> {
    // create A and B matrices for the Hamiltonian

    let mut A = Mat::<C64>::zeros(n_sites, n_sites);
    let mut B = Mat::<C64>::zeros(n_sites, n_sites);

    for c in couplings {
        let phase_factor = ((2. * J * PI) * (q.transpose() * c.inter_site_vector.as_ref())).exp();
        let (i, j) = (c.index1, c.index2);

        A[(i, j)] += (z[i].transpose() * c.matrix.as_ref() * z[j].conjugate()) * phase_factor;
        B[(i, j)] += (z[i].transpose() * c.matrix.as_ref() * z[j].as_ref()) * phase_factor;
    }

    A = _component_mul(&A, spin_coefficients);
    B = _component_mul(&B, spin_coefficients);

    // create Hamiltonian as a block matrix (the stack! macro creates a block matrix)
    let A_minus_C: Mat<C64> = A.clone() - C;
    let A_conj_minus_C: Mat<C64> = A.adjoint() - C;
    let hamiltonian: Mat<C64> = concat![[A_minus_C, B], [B.adjoint(), A_conj_minus_C]];

    // take square root of Hamiltonian using the LDL decomposition
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

    // calculate eigenvalues (energies) of the Hamiltonian and return
    (sqrt_hamiltonian.adjoint() * shc)
        .self_adjoint_eigenvalues(Side::Lower)
        .expect("Could not calculate eigenvalues of the Hamiltonian.")
}
/// Get the components of the rotation matrices for the axis indexed by `index`.
fn _get_rotation_component(rotations: &Vec<MatRef<C64>>, index: usize) -> Vec<Col<C64>> {
    rotations
        .par_iter()
        .map(|r| r.col(index).to_owned())
        .collect()
}

/// Perform componentwise multiplication on two matrices.
fn _component_mul(a: &Mat<C64>, b: &MatRef<C64>) -> Mat<C64> {
    let mut product = Mat::<C64>::zeros(a.nrows(), a.ncols());
    zip!(&mut product, a, b).for_each(|unzip!(product, x, y)| *product = x * y);
    product
}
