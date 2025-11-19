use crate::constants::SCALAR_J;
use crate::C64;
use faer::traits::ComplexField;
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

/// Create a block matrix from four sub-matrices.
#[inline]
pub fn block_matrix<T: ComplexField>(TL: &Mat<T>, TR: &Mat<T>, BL: &Mat<T>, BR: &Mat<T>) -> Mat<T> {
    let n_rows = TL.nrows() + BL.nrows();
    let n_cols = TL.ncols() + TR.ncols();
    let mut result = Mat::<T>::zeros(n_rows, n_cols);

    result
        .submatrix_mut(0, 0, TL.nrows(), TL.ncols())
        .copy_from(TL);
    result
        .submatrix_mut(0, TL.ncols(), TR.nrows(), TR.ncols())
        .copy_from(TR);
    result
        .submatrix_mut(TL.nrows(), 0, BL.nrows(), BL.ncols())
        .copy_from(BL);
    result
        .submatrix_mut(TL.nrows(), TR.ncols(), BR.nrows(), BR.ncols())
        .copy_from(BR);

    result
}
