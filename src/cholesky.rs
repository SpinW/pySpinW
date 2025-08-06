use nalgebra::{ComplexField, DMatrix};
use num_complex::Complex;

type C64 = Complex<f64>;

/// Cholesky decomposition for a complex matrix.
/// Modification of nalgebra's Cholesky::new() to fix a bug with complex matrices;
/// get rid of this when https://github.com/dimforge/nalgebra/issues/1536 is resolved
pub fn cholesky(mut matrix: DMatrix<C64>) -> Option<DMatrix<C64>> {
    assert!(matrix.is_square(), "The input matrix must be square.");         
    let n = matrix.nrows();         
    for j in 0..n {             
        for k in 0..j {                 
            let factor = unsafe { -matrix.get_unchecked((j, k)).clone() };                 
            let (mut col_j, col_k) = matrix.columns_range_pair_mut(j, k);                 
            let mut col_j = col_j.rows_range_mut(j..);                 
            let col_k = col_k.rows_range(j..);                 
            col_j.axpy(factor.conjugate(), &col_k, Complex::from(1.));             
        }
        // `diag` below is the diagonal elements of the L-matrix
        // if the input matrix is PD, then sqrt(diag) should always be real and positive
        // so we add a check to make sure this is the case
        let sqrt_denom = |v: C64| {
            if v == Complex::ZERO {                     
                return None;                 
            }                 
            match v.try_sqrt() {
                None => None,
                Some(sqrt) => {println!("{}", sqrt); if sqrt.im > f64::EPSILON {None} else {Some(sqrt)} }
            }
        };             
        let diag = unsafe { matrix.get_unchecked((j, j)).clone() };             
        if let Some(denom) = sqrt_denom(diag)            
        {                 
            unsafe {                     
                *matrix.get_unchecked_mut((j, j)) = denom.clone();                 
            }                 
            let mut col = matrix.view_range_mut(j + 1.., j);                 
            col /= denom;                 
            continue;             
        }             
        // The diagonal element is either zero or its square root could not             
        // be taken (e.g. for negative real numbers).             
        return None;         
    }
    Some(matrix)
}
