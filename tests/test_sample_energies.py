""" Samples give the magnon energies along a path """

import numpy as np

from pyspinw.exchange import HeisenbergExchange
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.path import Path
from pyspinw.sample import CrystalDomain, Multidomain, SingleCrystal
from pyspinw.site import Site
from pyspinw.structure import Structure
from pyspinw.symmetry.unitcell import UnitCell


def ferromagnetic_chain() -> Hamiltonian:
    """ Chain along x with J = -1, S = 1, so E(h) = 2 (1 - cos 2 pi h) meV, independent of k and l """
    site = Site(0, 0, 0, 0, 0, 1, name="S")
    return Hamiltonian(Structure([site], UnitCell(3.0, 8.0, 8.0)),
                       [HeisenbergExchange(site, site, j=-1.0, cell_offset=(1, 0, 0))])


PATH = Path([[0, 0, 0], [0.5, 0, 0]], n_points_per_segment=11)


def test_single_crystal_energies_are_the_hamiltonians():
    """ Check that a single crystal gives the Hamiltonian's own energies along the path """
    energies = SingleCrystal(ferromagnetic_chain()).energies(PATH)
    reference = ferromagnetic_chain().sorted_positive_energies(PATH)

    np.testing.assert_allclose(energies, reference)


def test_multidomain_energies_follow_each_domains_rotation():
    """ Check that each domain's transformation is applied to every q-point along the path.

    The unrotated domain gives the chain's dispersion along h. Rotated by 90 degrees about z, the
    path along h becomes one along k, across the chain, where the energy is flat.
    """
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])

    unrotated = Multidomain(ferromagnetic_chain(), [CrystalDomain(np.eye(3), 1.0)]).energies(PATH)
    rotated = Multidomain(ferromagnetic_chain(), [CrystalDomain(rotation, 1.0)]).energies(PATH)

    h = PATH.q_points()[:, 0]
    np.testing.assert_allclose(unrotated[0], 2 * (1 - np.cos(2 * np.pi * h)), atol=1e-9)
    np.testing.assert_allclose(rotated[0], 0.0, atol=1e-9)
