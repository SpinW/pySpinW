use nalgebra::{DMatrix, DVector};
use num_complex::Complex;
use std::cmp::max;

type C64 = Complex<f64>;

use lapack::zheev;
extern crate lapack_src;

/// Calculate eigenvalues and optionally eigenvectors for a square complex Hermitian matrix.
/// Wrapper around the LAPACK `zheev` routine.
/// See `zheev` docs here: https://netlib.org/lapack/explore-html/d8/d1c/group__heev_gadbb2b87ce42e51fdaac3228a857a58c8.html
pub fn eigs(
    mut matrix: DMatrix<C64>,
    vectors: bool,
) -> Result<nalgebra::DVector<f64>, &'static str> {
    // if jobz = 'V', vectors are calculated`
    let jobz = if vectors { b'V' } else { b'N' };

    // if uplo = 'U', the upper triangle of the matrix is used
    // (as matrix is Hermitian only one triangle is needed)
    let uplo = b'U';

    // n is the dimension of the matrix
    // for our purposes lda is the same
    let n = matrix.shape().0 as i32;
    // m is the matrix in lapack-friendly format
    let m: &mut [C64] = matrix.as_mut_slice();

    // w is the vector that the eigenvalues will be written to
    let mut output_vector = DVector::<f64>::zeros(n as usize);
    let w: &mut [f64] = output_vector.as_mut_slice();

    // unused here but required for the routine
    let mut info = 0;
    let mut rwork = vec![0.; max(1, 3 * (n as usize) - 2)];

    // if `lwork = -1`, `zheev` just calculates the optimal workspace size
    let mut placeholder = [Complex::from(0.)];
    unsafe {
        zheev(
            b'N',
            uplo,
            n,
            m,
            n,
            w,
            &mut placeholder,
            -1,
            &mut rwork,
            &mut info,
        )
    }

    match info {
        0 => (),
        x if x < 0 => return Err("LAPACK error: illegal argument."),
        x if x > 0 => return Err("LAPACK error: eigenvalue algorithm failed to converge"),
        _ => (),
    };

    let lwork = placeholder[0].re as i32;

    let mut workspace = vec![Complex::from(0.); lwork as usize];
    unsafe {
        zheev(
            jobz,            // whether to calculate eigenvectors
            uplo,            // whether to use upper or lower triangle of matrix
            n,               // size of matrix
            m,               // the matrix itself (gets overwritten with eigenvectors!!!)
            n,               // size of matrix (again)
            w,               // array eigenvalues are output into
            &mut workspace,  // workspace array for algorithm
            lwork,           // length of the workspace array
            &mut rwork,      // double precision array (don't know what this is for)
            &mut info,       // debug info output; <0 for invalid arguments, >0 for failed
                             //                    convergence, 0 for success
        )
    }

    match info {
        x if x < 0 => Err("LAPACK error: illegal argument."),
        x if x > 0 => Err("LAPACK error: eigenvalue algorithm failed to converge"),
        0 => Ok(output_vector),
        _ => unreachable!(),
    }
}
