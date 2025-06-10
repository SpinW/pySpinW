use std::f64::consts::PI;
use std::{iter::zip, ops::Sub};

use nalgebra::{stack, Const, DMatrix, DMatrixView, DVector, Dyn, Matrix3, Vector3};
use num_complex::Complex;
use numpy::{PyArray1, PyReadonlyArray1, PyReadonlyArray2, PyReadwriteArray2, ToPyArray};
use pyo3::prelude::*;
use rayon::iter::IntoParallelRefIterator;
use rayon::prelude::{IntoParallelIterator, ParallelIterator};

pub mod ldl;
use crate::ldl::ldl;

// some convenience types and statics for complex arithmetic
type C64 = Complex<f64>;
static J: C64 = Complex::new(0., 1.);

/// Temporary description of the coupling between atoms.
#[pyclass(frozen)]
pub struct Coupling {
    index1: usize,
    index2: usize,
    matrix: Matrix3<C64>,
    inter_site_vector: Vector3<f64>,
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
            matrix: matrix
                .try_as_matrix::<Const<3>, Const<3>, Dyn, Dyn>()
                .unwrap()
                .into(),
            inter_site_vector: inter_site_vector
                .try_as_matrix::<Const<3>, Const<1>, Dyn, Dyn>()
                .unwrap()
                .into(),
        }
    }
}

/// Run the main calculation step for a spinwave calculation.
#[pyfunction]
pub fn spinwave_calculation<'py>(
    py: Python<'py>,
    rotations: Vec<PyReadwriteArray2<C64>>,
    magnitudes: Vec<f64>,
    q_vectors: Vec<Vec<f64>>,
    couplings: Vec<Py<Coupling>>,
) -> PyResult<Vec<Bound<'py, PyArray1<C64>>>> {
    // convert PyO3-friendly array types to nalgebra matrices
    let r: Vec<Matrix3<C64>> = rotations
        .into_iter()
        .map(|m| -> Matrix3<C64> {
            let mv: DMatrixView<C64> = m.try_as_matrix().unwrap();
            mv.fixed_resize::<3, 3>(Complex::from(0.))
        })
        .collect();
    let qv = q_vectors.into_par_iter().map(Vector3::from_vec).collect();

    let c = couplings.par_iter().map(pyo3::Py::get).collect();

    let energies = _calc_spinwave(r, magnitudes, qv, c);
    Ok(energies.into_iter().map(|v| v.to_pyarray(py)).collect())
}

/// Run the main calculation step for a spinwave calculation.
#[allow(non_snake_case)]
fn _calc_spinwave(
    rotations: Vec<Matrix3<C64>>,
    magnitudes: Vec<f64>,
    q_vectors: Vec<Vector3<f64>>,
    couplings: Vec<&Coupling>,
) -> Vec<Vec<C64>> {
    let n_sites = rotations.len();

    // decompose rotation matrices
    // in the notation of Petit (2011)
    // eta[i] is the direction of the i'th moment in Cartesian coordinates
    let z: Vec<Vector3<C64>> = zip(
        _get_rotation_component(&rotations, 0),
        _get_rotation_component(&rotations, 1),
    )
    .map(|(r1, r2)| r1 + (r2 * J))
    .collect();
    let etas = _get_rotation_component(&rotations, 2);

    // make spin coefficients array
    // so spin_coefficients[i, j] = sqrt(S_i S_j) / 2
    let root_mags = DVector::<C64>::from_iterator(
        n_sites,
        magnitudes.iter().map(|x| Complex::from((0.5 * x).sqrt())),
    );
    let spin_coefficients = (root_mags.clone() * root_mags.transpose()).transpose();

    // create matrix C of Hamiltonian which is q-independent
    //
    // C_jj is the sum of eta_j^T * M * S_l eta_l over l
    // where M is the coupling matrix and S_i is the i'th spin
    // and we factor the l-dependent part out into `sites_term`
    // to avoid recalculating it for every coupling
    let sites_term: Vector3<C64> = zip(magnitudes, etas.clone())
        .map(|(magnitude, eta)| eta * Complex::from(magnitude))
        .sum::<Vector3<C64>>();
    let mut C = DMatrix::<C64>::zeros(n_sites, n_sites);
    for c in &couplings {
        *C.index_mut((c.index2, c.index2)) +=
            (etas[c.index2].transpose() * c.matrix * sites_term).into_scalar();
    }

    q_vectors
        .into_par_iter()
        .map(|q| _spinwave_single_q(q, &C, n_sites, &z, &spin_coefficients, &couplings))
        .collect()
}

/// Calculate spectra for a single q-value.
#[allow(non_snake_case)]
fn _spinwave_single_q(
    q: Vector3<f64>,
    C: &DMatrix<C64>,
    n_sites: usize,
    z: &[Vector3<C64>],
    spin_coefficients: &DMatrix<C64>,
    couplings: &Vec<&Coupling>,
) -> Vec<C64> {
    // create A and B matrices for the Hamiltonian

    let mut A = DMatrix::<C64>::zeros(n_sites, n_sites);
    let mut B = DMatrix::<C64>::zeros(n_sites, n_sites);

    for c in couplings {
        let phase_factor = ((2. * J * PI) * q.dot(&c.inter_site_vector)).exp();
        let (i, j) = (c.index1, c.index2);

        A[(i, j)] += (z[i].transpose() * c.matrix * z[j].conjugate()).into_scalar() * phase_factor;
        B[(i, j)] += (z[i].transpose() * c.matrix * z[j]).into_scalar() * phase_factor;
    }

    A = A.component_mul(spin_coefficients);
    B = B.component_mul(spin_coefficients);

    // create Hamiltonian as a block matrix (the stack! macro creates a block matrix)
    let A_minus_C: DMatrix<C64> = A.clone().sub(C);
    let A_conj_minus_C: DMatrix<C64> = A.adjoint().sub(C);
    let hamiltonian: DMatrix<C64> = stack![ A_minus_C, B; 
                                            B.adjoint(), A_conj_minus_C];

    // take square root of Hamiltonian using the LDL decomposition 
    let sqrt_hamiltonian = {
        let (l, d) = ldl(hamiltonian);
        l * DMatrix::from_diagonal(&d.map(nalgebra::Complex::sqrt))
    };

    // 'shc' is "square root of Hamiltonian with commutation"
    // We need to enforce the bosonic commutation properties, we do this
    // by finding the 'square root' of the matrix (i.e. finding K such that KK^dagger = H)
    // and then negating the second half.
    //
    // In matrix form we do
    //
    //     M = K^dagger g K
    //
    // where g is a diagonal matrix of length 2n, with the first n entries being 1, and the
    // remaining entries being -1. We do this by just multiplying the >n_sites columns of shc.
    let mut shc: DMatrix<C64> = sqrt_hamiltonian.clone();
    let mut negative_half = shc.view_mut((0, n_sites), (2 * n_sites, n_sites));
    negative_half *= Complex::from(-1.);

    // calculate eigenvalues (energies) of the Hamiltonian and return
    if let Some(v) = (shc.adjoint() * sqrt_hamiltonian).eigenvalues() {
        v.data.into()
    } else {
        panic!("Could not calculate eigenvalues of the Hamiltonian.")
    }
}
/// Get the components of the rotation matrices for the axis indexed by `index`.
fn _get_rotation_component(rotations: &Vec<Matrix3<C64>>, index: usize) -> Vec<Vector3<C64>> {
    rotations
        .par_iter()
        .map(|r| r.row(index).transpose())
        .collect()
}

/// A Python module implemented in Rust.
#[pymodule]
fn rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(spinwave_calculation, m)?)?;
    m.add_class::<Coupling>()?;
    Ok(())
}
