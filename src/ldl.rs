use nalgebra::{Const, DMatrix, DVector};
use num_complex::Complex;

type C64 = Complex<f64>;

// This algorithm is a modified copy of https://github.com/dimforge/nalgebra/pull/1515
// and when that pull request is merged, should be deleted and replaced with a call to
// `ldl` from the nalgebra library.

/// Computes the LDL^T factorization.
///
/// The input matrix `p` is assumed to be Hermitian and this decomposition will only read
/// the lower-triangular part of `p`.
pub fn ldl(p: DMatrix<C64>) -> (DMatrix<C64>, DVector<C64>) {
    let n = p.ncols();

    let n_dim = p.shape_generic().1;

    let mut d = DVector::<C64>::zeros_generic(n_dim, Const::<1>);
    let mut l = DMatrix::<C64>::zeros_generic(n_dim, n_dim);

    for j in 0..n {
        let mut d_j = p[(j, j)];

        if j > 0 {
            for k in 0..j {
                d_j -= d[k] * l[(j, k)].clone().powi(2);
            }
        }

        d[j] = d_j;

        for i in j..n {
            let mut l_ij = p[(j, i)];

            for k in 0..j {
                l_ij -= d[k] * l[(j, k)] * l[(i, k)];
            }

            if d[j] == Complex::from(0.) {
                l[(i, j)] = Complex::from(0.)
            } else {
                l[(i, j)] = l_ij / d[j];
            }
        }
    }

    (l, d)
}
