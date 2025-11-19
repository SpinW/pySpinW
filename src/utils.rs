use crate::constants::SCALAR_J;
use crate::C64;
use faer::{unzip, zip, Col, ColRef, Mat, MatRef};
/// Utility functions for the calculations.
use rayon::iter::{IntoParallelIterator, ParallelIterator};

/// Get the components of the rotation matrices for the axis indexed by `index`.
#[inline(always)]
pub fn get_rotation_components(rotations: Vec<MatRef<C64>>) -> (Vec<Col<C64>>, Vec<ColRef<C64>>) {
    // r.col(index) gets the components of the rotation matrix
    // and then in the map we compile the x and y components into z, and the other into eta
    // so this function returns (z, eta)
    rotations
        .into_par_iter()
        .map(|r| (r.col(0) + (r.col(1) * SCALAR_J), r.col(2)))
        .collect()
}

/// Perform componentwise multiplication on two matrices.
#[inline]
pub fn component_mul(a: &Mat<C64>, b: &Mat<C64>) -> Mat<C64> {
    let mut product = Mat::<C64>::zeros(a.nrows(), a.ncols());
    zip!(&mut product, a, b).for_each(|unzip!(product, x, y)| *product = x * y);
    product
}

/// Combine the A, B, and C matrices into the Hamiltonian `h(q)` matrix.
#[inline(always)]
pub fn make_block_hamiltonian(A: Mat<C64>, B: Mat<C64>, C: &Mat<C64>) -> Mat<C64> {
    let A_minus_C: Mat<C64> = A.clone() - C;
    let A_conj_minus_C: Mat<C64> = A.adjoint() - C;

    let n_sites = B.nrows();

    let mut hamiltonian = Mat::<C64>::zeros(2 * n_sites, 2 * n_sites);

    // `submatrix_mut` takes a writable view of part of the matrix, and `copy_from` copies
    // its argument into the matrix. So this is copying each of our matrices into sections
    // of the `hamiltonian` matrix
    // copy A - C into upper-left quarter of matrix
    hamiltonian
        .submatrix_mut(0, 0, n_sites, n_sites)
        .copy_from(A_minus_C);
    // copy B into upper-right
    hamiltonian
        .submatrix_mut(n_sites, 0, n_sites, n_sites)
        .copy_from(B.as_ref());
    // copy B* into bottom-left
    hamiltonian
        .submatrix_mut(0, n_sites, n_sites, n_sites)
        .copy_from(B.adjoint());
    // copy A* - C into bottom-right
    hamiltonian
        .submatrix_mut(n_sites, n_sites, n_sites, n_sites)
        .copy_from(A_conj_minus_C);

    hamiltonian
}
