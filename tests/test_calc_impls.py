"""Compare spinwave calculation implementations.

This file checks that the the Python and Rust implementations of
the spinwave calculation both return the same results for a set
of examples.

"""

import numpy as np
import pytest

from examples.raw_calculations.utils import py_classes, InternalClasses

try:
    from pyspinw.rust import spinwave_calculation as rs_spinwave, Coupling as RsCoupling, MagneticField as RsField
    rs_classes = InternalClasses(coupling=RsCoupling, field=RsField)
except ImportError:
    # we can use the --runxfail option for pytest to then ensure
    # that the Rust tests run and pass if we're expecting Rust to be installed
    pytestmark = pytest.mark.xfail(raises=NameError, reason="Rust module not installed.")

from pyspinw.calculations.spinwave import spinwave_calculation as py_spinwave

from examples.raw_calculations.ferromagnetic_chain import heisenberg_ferromagnet
from examples.raw_calculations.ferromagnet_gtensor import ferromagnet_gtensor
from examples.raw_calculations.antiferro_chain import antiferro_chain
from examples.raw_calculations.antiferro_ef import antiferro_ef
from examples.raw_calculations.kagome import kagome_ferromagnet
from examples.raw_calculations.kagome_antiferro import kagome_antiferromagnet
from examples.raw_calculations.kagome_supercell import kagome_supercell

@pytest.mark.rust
@pytest.mark.parametrize("example",
                         [heisenberg_ferromagnet,
                          ferromagnet_gtensor,
                          antiferro_chain,
                          antiferro_ef,
                          kagome_ferromagnet,
                          kagome_antiferromagnet,
#                          kagome_supercell,
                          ])
def test_calc_impls(example):
    """Compare Rust and Python spinwave calculation implementations."""
    rs_energies, rs_sqw = rs_spinwave(*example(classes=rs_classes))
    py_energies, py_sqw = py_spinwave(*example(classes=py_classes))

    # we test to an absolute tolerance of 1e-6 in line with the MATLAB
    np.testing.assert_allclose(np.sort(rs_energies), np.sort(py_energies), atol=1e-6, rtol=0)
    np.testing.assert_allclose(np.sort(rs_sqw), np.sort(py_sqw), atol=1e-6, rtol=1e-3)
