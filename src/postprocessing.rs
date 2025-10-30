/// Post-processing functions for spin wave calculations.
use faer::Mat;

use crate::C64;

/// Calculate the neutron scattering cross-section.
///
/// # Parameters
/// `Sab`: The correlation tensor returned by `calc_spinwave`.
/// This is a vector over q, where each element is a vector over (non-zero) omega where each
/// element is a 3x3 matrix of complex numbers representing S^alpha,beta(q, omega).
pub fn neutron(Sab: Vec<Vec<Mat<C64>>>) -> Vec<Vec<Mat<C64>>> {
    todo!()
}
