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
use crate::spinwave::{calc_energies, calc_spinwave};

mod constants;
mod utils;

// convenience type for complex arithmetic
type C64 = Complex<f64>;

// nicer names for PyO3 output types
type Energies<'py> = Vec<Bound<'py, PyArray1<f64>>>;
type SabTensor<'py> = Vec<Vec<Bound<'py, PyArray2<C64>>>>;
type SQw<'py> = Vec<Bound<'py, PyArray1<f64>>>;

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

/// Calculate the energies (eigenvalues) for a system.
///
/// # Parameters
/// - `rotations`: A list of 2D numpy arrays representing rotation matrices for each atom.
/// - `magnitudes`: A list of magnitudes for each atom.
/// - `q_vectors`: A list of q-vectors where the energies should be calculated.
/// - `couplings`: A list of `Coupling` objects representing the interactions between atoms
/// - `field`: An optional `MagneticField` object representing an external magnetic field.
///
/// # Returns
/// A list of 1D numpy arrays, each containing the energies for the corresponding q-vector.
#[pyfunction(signature = (rotations, magnitudes, q_vectors, couplings, field=None))]
pub fn energies<'py>(
    py: Python<'py>,
    rotations: Vec<PyReadonlyArray2<C64>>,
    magnitudes: Vec<f64>,
    q_vectors: Vec<Vec<f64>>,
    couplings: Vec<Py<Coupling>>,
    field: Option<MagneticField>,
) -> PyResult<Energies<'py>> {
    // convert PyO3-friendly array types to faer matrices
    let r: Vec<MatRef<C64>> = rotations
        .into_iter()
        .map(faer_ext::IntoFaer::into_faer)
        .collect();

    let c = couplings.par_iter().map(pyo3::Py::get).collect();

    let results = calc_energies(r, magnitudes, q_vectors, c, field);
    Ok(results
        .into_iter()
        .map(|result| result.to_pyarray(py))
        .collect())
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
/// - `field`: An optional `MagneticField` object representing an external magnetic field.
///
/// # Returns
/// A tuple containing:
/// - A list of 1D numpy arrays, each containing the energies for the corresponding q-vector.
/// - A list of 1D numpy arrays, each containing the neutron scattering cross-section
///   for the corresponding q-vector (indexed over omega).
#[pyfunction(signature = (rotations, magnitudes, q_vectors, couplings, positions, field=None))]
pub fn spinwave_calculation<'py>(
    py: Python<'py>,
    rotations: Vec<PyReadonlyArray2<C64>>,
    magnitudes: Vec<f64>,
    q_vectors: Vec<Vec<f64>>,
    couplings: Vec<Py<Coupling>>,
    positions: Vec<PyReadonlyArray1<f64>>,
    field: Option<MagneticField>,
) -> PyResult<(Energies<'py>, SQw<'py>)> {
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

    let results = calc_spinwave(r, magnitudes, q_vectors.clone(), c, p, field, false);
    Ok((
        results
            .iter()
            .map(|result| result.energies.to_pyarray(py))
            .collect(),
        results
            .iter()
            .map(|result| result.intensities.to_pyarray(py))
            .collect(),
    ))
}

/// Same as spinwave_calculation but also returns Sab tensors.
#[pyfunction(signature = (rotations, magnitudes, q_vectors, couplings, positions, field=None))]
pub fn spinwave_calculation_Sab<'py>(
    py: Python<'py>,
    rotations: Vec<PyReadonlyArray2<C64>>,
    magnitudes: Vec<f64>,
    q_vectors: Vec<Vec<f64>>,
    couplings: Vec<Py<Coupling>>,
    positions: Vec<PyReadonlyArray1<f64>>,
    field: Option<MagneticField>,
) -> PyResult<(Energies<'py>, SQw<'py>, SabTensor<'py>)> {
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

    let results = calc_spinwave(r, magnitudes, q_vectors.clone(), c, p, field, true);
    Ok((
        results
            .iter()
            .map(|result| result.energies.to_pyarray(py))
            .collect(),
        results
            .iter()
            .map(|result| result.intensities.to_pyarray(py))
            .collect(),
        results
            .iter()
            .map(|result| {
                let sab = result.sab.clone().unwrap();
                sab.iter()
                    .map(|sab_qw| PyArray2::from_array(py, &sab_qw.as_ref().into_ndarray()))
                    .collect()
            })
            .collect(),
    ))
}

/// A Python module implemented in Rust.
#[pymodule]
fn rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(energies, m)?)?;
    m.add_function(wrap_pyfunction!(spinwave_calculation, m)?)?;
    m.add_function(wrap_pyfunction!(spinwave_calculation_Sab, m)?)?;
    m.add_class::<Coupling>()?;
    m.add_class::<MagneticField>()?;
    Ok(())
}
