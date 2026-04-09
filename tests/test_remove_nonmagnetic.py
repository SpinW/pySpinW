""" Tests for the removal of nonmagnetic sites """
import numpy as np
from pyspinw import UnitCell, LatticeSite, Structure, HeisenbergExchange, AxisMagnitudeAnisotropy, Hamiltonian


def test_without_nonmagnetic_structure():
    """ Check that removing nonmagnetic sites from structures works """
    unit_cell = UnitCell(1,1,1)

    magnetic = LatticeSite(0.1,0,0,0,0,1)
    nonmagnetic = LatticeSite(0.2,0,0,0,0,0)

    structure = Structure([magnetic, nonmagnetic], unit_cell=unit_cell)

    without = structure.without_nonmagnetic()

    for site in without.sites:
        assert not np.allclose(site.spin_data, 0.0)


def test_without_nonmagnetic_hamiltonian():
    """ Check that removing nonmagnetic sites from hamiltonian works

    Assumes that structure removal works (tested elsewhere)
    """

    unit_cell = UnitCell(1, 1, 1)

    magnetic = LatticeSite(0.1, 0, 0, 0, 0, 1)
    nonmagnetic = LatticeSite(0.2, 0, 0, 0, 0, 0)

    exchanges = [
        HeisenbergExchange(magnetic, magnetic, j=1),
        HeisenbergExchange(magnetic, nonmagnetic, j=1),
        HeisenbergExchange(nonmagnetic, nonmagnetic, j=1),
        HeisenbergExchange(nonmagnetic, magnetic, j=1)]

    anisotropies = [
        AxisMagnitudeAnisotropy(magnetic, 1, [1,0,0]),
        AxisMagnitudeAnisotropy(nonmagnetic, 1, [1,0,0])
    ]

    structure = Structure([magnetic, nonmagnetic], unit_cell=unit_cell)

    hamiltonian = Hamiltonian(structure, exchanges, anisotropies)

    without = hamiltonian.without_nonmagnetic()

    for exchange in without.exchanges:
        assert not np.allclose(exchange.site_1.spin_data, 0.0)
        assert not np.allclose(exchange.site_2.spin_data, 0.0)

    for anisotropy in without.anisotropies:
        assert not np.allclose(anisotropy.site.spin_data, 0.0)
