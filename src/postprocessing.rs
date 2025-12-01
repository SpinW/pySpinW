/// Post-processing functions for spin wave calculations.
use faer::{Col, Mat};
use rayon::prelude::*;

use crate::utils::component_mul;
use crate::C64;

/// Calculate the neutron scattering cross-section.
///
/// # Parameters
/// `Sab`: The correlation tensor returned by `calc_spinwave`.
/// This is a vector over q, where each element is a vector over (non-zero) omega where each
/// element is a 3x3 matrix of complex numbers representing S^alpha,beta(q, omega).
pub fn neutron(Sab: Vec<Vec<Mat<C64>>>, q_vectors: Vec<Vec<f64>>) -> Vec<Vec<f64>> {
    Sab.par_iter()
        .zip(q_vectors)
        .map(|(Sab_q, q)| {
            Sperp_single_q(
                &Sab_q,
                Col::<C64>::from_iter(q.iter().map(|x| C64::from(x))),
            )
        })
        .collect()
}

fn Sperp_single_q(Sab_q: &Vec<Mat<C64>>, wavevector: Col<C64>) -> Vec<f64> {
    let mut norm_q = wavevector.as_ref() / wavevector.norm_l2();
    if norm_q.has_nan() {
        norm_q = Col::<C64>::from_iter(vec![C64::from(0.0), C64::from(0.0), C64::from(0.0)]);
    }
    let perp_factor = Mat::<C64>::identity(3, 3) - (norm_q.as_ref() * norm_q.adjoint());
    Sab_q
        .iter()
        .map(|Sab_qw| -> f64 { component_mul(Sab_qw, &perp_factor).sum().re })
        .collect()
}
