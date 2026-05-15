//! A Rust-based implementation of the spinwave calculation.
//! The file `lib.rs` contains the Python bindings,
//! and the actual pure Rust calculation is in `spinwave.rs`.
#![allow(non_snake_case)]

use faer::{Col, ColRef, Mat, MatRef};
use faer_ext::{IntoFaer, IntoNdarray};
use num_complex::Complex;
use numpy::{PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2, ToPyArray};
use pyo3::prelude::*;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

mod spinwave;
use crate::spinwave::calc_spinwave;

mod constants;
mod utils;

// convenience type for complex arithmetic
type C64 = Complex<f64>;

// nicer names for PyO3 output types
type Energies<'py> = Vec<Bound<'py, PyArray1<C64>>>;
type SabTensor<'py> = Vec<Vec<Bound<'py, PyArray2<C64>>>>;
type SQw<'py> = Vec<Bound<'py, PyArray1<f64>>>;

type Wavefunction<'py> = Vec<Bound<'py, PyArray2<C64>>>;

/// Temporary description of the coupling between atoms.
#[pyclass(frozen)]
pub struct Coupling {
    index1: usize,
    index2: usize,
    matrix: Mat<C64>,
    inter_site_vector: Col<f64>,
}

#[pymethods]
impl Coupling {
    #[new]
    fn new(
        index1: usize,
        index2: usize,
        matrix: PyReadonlyArray2<C64>,
        inter_site_vector: PyReadonlyArray1<f64>,
    ) -> Self {
        Coupling {
            index1,
            index2,
            matrix: matrix.into_faer().to_owned(),
            inter_site_vector: inter_site_vector.into_faer().to_owned(),
        }
    }
    fn __eq__(&self, other: &Self) -> bool {
        self == other
    }
}

impl PartialEq for Coupling {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.index1 == other.index1
            && self.index2 == other.index2
            && (self.matrix.clone() - other.matrix.clone()).norm_l1() < 1e-6
            && (self.inter_site_vector.clone() - other.inter_site_vector.clone()).norm_l1() < 1e-6
    }
}

#[pyclass(frozen)]
#[derive(Clone)]
pub struct MagneticField {
    vector: Col<C64>,
    g_tensors: Vec<Mat<C64>>,
}

#[pymethods]
impl MagneticField {
    #[new]
    fn new(vector: PyReadonlyArray1<C64>, g_tensors: Vec<PyReadonlyArray2<C64>>) -> Self {
        MagneticField {
            vector: vector.into_faer().to_owned(),
            g_tensors: g_tensors
                .into_iter()
                .map(|t| t.into_faer().to_owned())
                .collect(),
        }
    }
}

/// Calculate energies and neutron scattering cross-section for a system.
///
/// # Parameters
/// - `rotations`: A list of 2D numpy arrays representing rotation matrices for each atom.
/// - `magnitudes`: A list of magnitudes for each atom.
/// - `q_vectors`: A list of q-vectors where the calculations should be performed.
/// - `couplings`: A list of `Coupling` objects representing the interactions between atoms.
/// - `positions`: A list of 1D numpy arrays representing the relative positions of each atom
///  in the unit cell.
/// - `rlu_to_cart`: An optional 3x3 matrix to convert from reciprocal lattice units to Cartesian.
/// - `field`: An optional `MagneticField` object representing an external magnetic field.
///
/// # Returns
/// A tuple containing:
/// - A list of 1D numpy arrays, each containing the energies for the corresponding q-vector.
/// - A list of 1D numpy arrays, each containing the neutron scattering cross-section
///   for the corresponding q-vector (indexed over omega).
// #[pyfunction(signature = (rotations, magnitudes, q_vectors, couplings, positions, rlu_to_cart=None, field=None, rotating_frame=None))]
// pub fn spinwave_calculation<'py>(
//     py: Python<'py>,
//     rotations: Vec<PyReadonlyArray2<C64>>,
//     magnitudes: Vec<f64>,
//     q_vectors: Vec<Vec<f64>>,
//     couplings: Vec<Py<Coupling>>,
//     positions: Vec<PyReadonlyArray1<f64>>,
//     rlu_to_cart: Option<PyReadonlyArray2<f64>>,
//     field: Option<MagneticField>,
//     rotating_frame: Option<Vec<PyReadonlyArray1<f64>>>,
// ) -> PyResult<(Energies<'py>, SQw<'py>)> {
//     // convert PyO3-friendly array types to faer matrices
//     let r: Vec<MatRef<C64>> = rotations
//         .into_iter()
//         .map(faer_ext::IntoFaer::into_faer)
//         .collect();
//
//     let c = couplings.par_iter().map(pyo3::Py::get).collect();
//
//     let p: Vec<ColRef<f64>> = positions
//         .into_iter()
//         .map(faer_ext::IntoFaer::into_faer)
//         .collect();
//
//     let to_cart: Option<MatRef<f64>> = rlu_to_cart.map(|f| f.into_faer());
//     let rotating: Option<Vec<ColRef<f64>>> =
//         rotating_frame.map(|f| f.into_iter().map(faer_ext::IntoFaer::into_faer).collect());
//
//     let results = calc_spinwave(
//         r,
//         magnitudes,
//         q_vectors.clone(),
//         c,
//         p,
//         to_cart,
//         field,
//         rotating,
//         false,
//     );
//     Ok((
//         results
//             .iter()
//             .map(|result| result.energies.to_pyarray(py))
//             .collect(),
//         results
//             .iter()
//             .map(|result| result.intensities.to_pyarray(py))
//             .collect(),
//     ))
// }

/// Same as spinwave_calculation but also returns Sab tensors.
#[pyfunction(signature = (rotations, magnitudes, q_vectors, couplings, positions, rlu_to_cart=None, field=None, rotating_frame=None, save_sab=false, save_wavefunctions=false))]
pub fn spinwave_calculation<'py>(
    py: Python<'py>,
    rotations: Vec<PyReadonlyArray2<C64>>,
    magnitudes: Vec<f64>,
    q_vectors: Vec<Vec<f64>>,
    couplings: Vec<Py<Coupling>>,
    positions: Vec<PyReadonlyArray1<f64>>,
    rlu_to_cart: Option<PyReadonlyArray2<f64>>,
    field: Option<MagneticField>,
    rotating_frame: Option<Vec<PyReadonlyArray1<f64>>>,
    save_sab: bool,
    save_wavefunctions: bool,
) -> PyResult<(Energies<'py>, SQw<'py>, Option<SabTensor<'py>>, Option<Wavefunction<'py>>)> {
    // convert PyO3-friendly array types to faer matrices
    let r: Vec<MatRef<C64>> = rotations
        .into_iter()
        .map(faer_ext::IntoFaer::into_faer)
        .collect();

    let c = couplings.par_iter().map(pyo3::Py::get).collect();

    let p: Vec<ColRef<f64>> = positions
        .into_iter()
        .map(faer_ext::IntoFaer::into_faer)
        .collect();

    let to_cart: Option<MatRef<f64>> = rlu_to_cart.map(|f| f.into_faer());
    let rotating: Option<Vec<ColRef<f64>>> =
        rotating_frame.map(|f| f.into_iter().map(faer_ext::IntoFaer::into_faer).collect());

    let results = calc_spinwave(
        r,
        magnitudes,
        q_vectors.clone(),
        c,
        p,
        to_cart,
        field,
        rotating,
        save_sab,
        save_wavefunctions,
    );
    Ok((
        // Energies, vector of arrays
        results
            .iter()
            .map(|result| result.energies.to_pyarray(py))
            .collect(),

        // Intensity, vector of arrays
        results
            .iter()
            .map(|result| result.intensities.to_pyarray(py))
            .collect(),

        // Sab, optional vector of vector of complex matrices
        match save_sab {
            false => None,
            true => Some(
                results
                    .iter()         // iterating of vector spinwave results
                    .map(|result| {      // result should be of type SpinwaveResult
                        result.sab.as_ref().map(|sab| {     // map over the optional sab entries, vectors of matrices
                            sab.iter()             // Convert to pythong objects
                                .map(|sab_qw| PyArray2::from_array(py, &sab_qw.as_ref().into_ndarray()))
                                .collect()
                        })
                    })
                    .flatten() // Remove Options
                    .collect()
                )
            },

        // Wavefunction, optional list of matrices
        match save_wavefunctions {
            false => None,
            true => Some(
                results
                    .iter()
                    .map(|result| {
                        result.wavefunctions.as_ref().map(|wavefunctions| {
                            PyArray2::from_array(py, &wavefunctions.as_ref().into_ndarray())
                        })
                    })
                    .flatten()
                    .collect(),
            )
        }
    ))
}

/// A Python module implemented in Rust.
#[pymodule]
fn rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(spinwave_calculation, m)?)?;
    //m.add_function(wrap_pyfunction!(spinwave_calculation_Sab, m)?)?;
    m.add_class::<Coupling>()?;
    m.add_class::<MagneticField>()?;
    Ok(())
}
