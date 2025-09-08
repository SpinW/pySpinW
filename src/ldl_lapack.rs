use nalgebra::{DMatrix, DVector};
use num_complex::Complex;

type C64 = Complex<f64>;

use lapack::zhetrf;
extern crate lapack_src;

/// Perform the LDL (Bunch-Kaufman) decomposition for a complex square Hermitian matrix.
/// Wrapper around the LAPACK `zhetrf` routine.
/// See `zhetrf` docs here: https://netlib.org/lapack/explore-html/d8/d0e/group__hetrf_ga114675389727c72322841cef488c10dd.html
pub fn ldl(mut matrix: DMatrix<C64>) -> (DMatrix<C64>, DVector<C64>) {
    // n is the dimension of the matrix
    // for our purposes lda is the same
    let n = matrix.shape().0 as i32;
    // m is the matrix in lapack-friendly format
    let m: &mut [C64] = matrix.as_mut_slice();

    // matrix giving pivot information (i.e. whether D contains 2x2 blocks)
    let mut ipiv = vec![0; n as usize];

    // we don't bother calculating LWORK - scipy seems to also just use N
    let mut workspace = vec![Complex::from(0.); n as usize];

    let mut info = 0;

    unsafe {
        zhetrf(
            b'L',            // whether to do LDL ('L') or UDU ('U')
            n,               // size of matrix (we assume it's square - this isn't necessary)
            m,               // the matrix itself (overwritten with decomposition!)
            n,               // size of matrix again
            &mut ipiv,       // info on pivots (i.e. whether matrix contains 2x2 blocks)
            &mut workspace,  // workspace array for algorithm
            n,               // size of workspace
            &mut info        // debug info output
            )
    };

    if info < 0 {panic!("LDL factorisation failed. (code {})", info)};

    // LAPACK returns the decomposition as a matrix which needs unpacking;
    // we just copy how SciPy does it in `scipy.linalg._decomp_ldl._ldl_get_d_and_l`
    // https://github.com/scipy/scipy/blob/v1.16.1/scipy/linalg/_decomp_ldl.py#L246
    // TODO: the above ^

    let d = matrix.diagonal();
    let l = matrix.lower_triangle();

    (l, d)
} 
