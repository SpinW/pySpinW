//! A Rust-based implementation of the spinwave calculation.
//! The file `lib.rs` contains the Python bindings,
//! and the actual pure Rust calculation is in `spinwave.rs`.

use faer::{Col, Mat, MatRef};
use faer_ext::IntoFaer;
use num_complex::Complex;
use numpy::{PyArray1, PyReadonlyArray1, PyReadonlyArray2, ToPyArray};
use pyo3::prelude::*;
use rayon::iter::{IntoParallelIterator, IntoParallelRefIterator, ParallelIterator};

pub mod spinwave;
use crate::spinwave::calc_spinwave;

// convenience type for complex arithmetic
type C64 = Complex<f64>;

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

/// Run the main calculation step for a spinwave calculation.
#[pyfunction]
pub fn spinwave_calculation<'py>(
    py: Python<'py>,
    rotations: Vec<PyReadonlyArray2<C64>>,
    magnitudes: Vec<f64>,
    q_vectors: Vec<Vec<f64>>,
    couplings: Vec<Py<Coupling>>,
) -> PyResult<Vec<Bound<'py, PyArray1<f64>>>> {
    // convert PyO3-friendly array types to faer matrices
    let r: Vec<MatRef<C64>> = rotations.into_iter().map(faer_ext::IntoFaer::into_faer).collect();
    let qv = q_vectors.into_par_iter().map(Col::from_iter).collect();

    let c = couplings.par_iter().map(pyo3::Py::get).collect();

    let energies = calc_spinwave(r, magnitudes, qv, c);
    Ok(energies.into_iter().map(|v| v.to_pyarray(py)).collect())
}

/// A Python module implemented in Rust.
#[pymodule]
fn rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(spinwave_calculation, m)?)?;
    m.add_class::<Coupling>()?;
    Ok(())
}
