""" A powder spectrum can be taken over an explicit list of |q| (issue #254) """

import numpy as np

from pyspinw.exchange import HeisenbergExchange
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.path import EmpiricalPath1D
from pyspinw.sample import Powder
from pyspinw.site import LatticeSite
from pyspinw.structure import Structure
from pyspinw.symmetry.unitcell import UnitCell


def test_empirical_path_counts_its_points():
    """ The path knows how many q values it holds """
    assert EmpiricalPath1D([0.37, 0.5, 0.8]).n_points == 3


def test_powder_spectrum_takes_a_plain_list_of_q():
    """ The signature accepts an ArrayLike, which is wrapped in an EmpiricalPath1D """
    site = LatticeSite(0, 0, 0, 0, 0, 1, name="S")
    hamiltonian = Hamiltonian(Structure([site], UnitCell(3.0, 8.0, 8.0)),
                              [HeisenbergExchange(site, site, j=-1.0, cell_offset=(1, 0, 0))])

    spectrum = Powder(hamiltonian).spectrum([0.37, 0.5], n_samples=100, n_energy_bins=10)

    assert np.array(spectrum[-1]).shape == (2, 10)
