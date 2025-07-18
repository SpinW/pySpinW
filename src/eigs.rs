use nalgebra::{DMatrix, DVector};
use num_complex::Complex;

type C64 = Complex<f64>;

use lapack::zheev;


/// Calculate eigenvalues and optionally eigenvectors for a square complex Hermitian matrix.
/// See `zheev` docs here: https://netlib.org/lapack/explore-html/d8/d1c/group__heev_gadbb2b87ce42e51fdaac3228a857a58c8.html
pub fn eigs(mut matrix: DMatrix<C64>, vectors: bool) -> Result<nalgebra::DVector<f64>, &'static str> {

    // if jobz = 'V', vectors are calculated`
    let jobz = if vectors {b'V'} else {b'N'};

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
    let w = output_vector.as_mut_slice();

    // unused here but required for the routine
    let mut info = 0;
    let mut placeholder = [0.];

    // if `lwork = -1`, `zheev` just calculates the optimal workspace size
    let mut workspace = [Complex::from(0.)];
    unsafe { zheev(b'N', uplo, n, m, n, w, &mut workspace, -1, &mut placeholder, &mut info) }

    match info {
        0 => (),
        x if x < 0 => return Err("LAPACK error: illegal argument."),
        x if x > 0 => return Err("LAPACK error: eigenvalue algorithm failed to converge"),
        _ => ()
    };
    
    let lwork = workspace[0].re as i32;

    let mut workspace = vec![Complex::from(0.); lwork as usize];
    unsafe { zheev(jobz, uplo, n, m, n, w, &mut workspace, lwork, &mut placeholder, &mut info) }

    match info {
        x if x < 0 => Err("LAPACK error: illegal argument."),
        x if x > 0 => Err("LAPACK error: eigenvalue algorithm failed to converge"),
        _ => Ok(output_vector)
    }
}
