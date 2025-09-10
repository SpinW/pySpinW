"""Compare spinwave calculation implementations.

This file checks that the the Python and Rust implementations of
the spinwave calculation both return the same results for a set
of examples.

"""

import numpy as np
import pytest

try:
    from pyspinw.rust import spinwave_calculation as rs_spinwave, Coupling as RsCoupling
except ImportError:
    pytestmark = pytest.mark.skip("Rust module not installed.")

from pyspinw.calculations.spinwave import spinwave_calculation as py_spinwave, Coupling as PyCoupling

from examples.raw_calculations.ferromagnetic_chain import heisenberg_ferromagnet
from examples.raw_calculations.antiferro_chain import antiferro_chain
from examples.raw_calculations.kagome import kagome_ferromagnet
from examples.raw_calculations.kagome_antiferro import kagome_antiferromagnet
from examples.raw_calculations.kagome_supercell import kagome_supercell

@pytest.mark.parametrize("example",
                         [heisenberg_ferromagnet,
                          antiferro_chain,
                          kagome_ferromagnet,
                          kagome_antiferromagnet,
                          kagome_supercell,])
def test_calc_impls(example):
    """Compare Rust and Python spinwave calculation implementations."""
    rs_results = rs_spinwave(*example(coupling_class=RsCoupling))
    py_results = py_spinwave(*example(coupling_class=PyCoupling))

    # we test to an absolute tolerance of 1e-6 in line with the MATLAB
    np.testing.assert_allclose(np.sort(rs_results), np.sort(py_results), atol=1e-6, rtol=0)
