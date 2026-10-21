""" Check that the classical minimiser reads the magnetic field in tesla, like the spin wave calculation """

import numpy as np

from pyspinw.calculations.energy_minimisation import ClassicalEnergyMinimisation
from pyspinw.constants import MU_B
from pyspinw.exchange import HeisenbergExchange
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.site import LatticeSite
from pyspinw.structure import Structure
from pyspinw.symmetry.unitcell import UnitCell


def _ferromagnetic_chain(sz: float) -> Hamiltonian:
    """ One-site ferromagnetic chain with S=1 along z, with sign sz """
    site = LatticeSite(0, 0, 0, 0, 0, sz, name="S")
    exchange = HeisenbergExchange(site, site, j=-1.0, cell_offset=(1, 0, 0))
    return Hamiltonian(Structure([site], UnitCell(3.0, 8.0, 8.0)), [exchange])


def _antiferromagnetic_chain() -> Hamiltonian:
    """ Two-site antiferromagnetic chain, J=1 meV, S=1, spins along z """
    a = LatticeSite(0, 0, 0, 0, 0, 1, name="A")
    b = LatticeSite(0.5, 0, 0, 0, 0, -1, name="B")
    exchanges = [HeisenbergExchange(a, b, j=1.0, cell_offset=(0, 0, 0)),
                 HeisenbergExchange(b, a, j=1.0, cell_offset=(1, 0, 0))]
    return Hamiltonian(Structure([a, b], UnitCell(6.0, 8.0, 8.0)), exchanges)


def test_zeeman_energy_is_in_tesla():
    """ E(up) - E(down) in a field B along z should be 2 g mu_B S B, with B in tesla """
    field = np.array([0.0, 0.0, 1.0])

    up = ClassicalEnergyMinimisation(_ferromagnetic_chain(1.0), field=field).energy()
    down = ClassicalEnergyMinimisation(_ferromagnetic_chain(-1.0), field=field).energy()

    g, spin, b = 2.0, 1.0, 1.0
    assert np.isclose(up - down, 2 * g * MU_B * spin * b)


def test_antiferromagnet_cants_by_the_tesla_angle():
    """ A transverse field B cants an antiferromagnetic chain by sin(theta) = g mu_B B / 4 J S, B in tesla

    This goes through the forces rather than the energy alone. The net moment along the field is what the field
    fixes; the staggered direction is soft and need not converge as tightly.
    """
    field = np.array([1.0, 0.0, 0.0])

    found = _antiferromagnetic_chain().ground_state(field=field, verbose=False, seed=1)

    g, spin, j, b = 2.0, 1.0, 1.0, 1.0
    sin_theta = g * MU_B * b / (4 * j * spin)
    net_moment_x = np.mean([site.base_spin[0] for site in found.structure.sites])

    # The spins tilt away from the field
    assert np.isclose(net_moment_x, -spin * sin_theta, atol=1e-5)
