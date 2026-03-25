import numpy as np
import pytest

from pyspinw.calculations.energy_minimisation import ClassicalEnergyMinimisation
from pyspinw.cell_offsets import CellOffset
from pyspinw.interface import generate_exchanges
from pyspinw.coupling import HeisenbergCoupling
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.site import LatticeSite
from pyspinw.structures import Structure
from pyspinw.symmetry.supercell import TrivialSupercell, SummationSupercell, CommensuratePropagationVector, \
    TransformationSupercell, RotationTransform

from pyspinw.symmetry.unitcell import UnitCell

def test_energy_invariance_trivial_supercell():
    """ Check energy stays the same when we make the supercell bigger """
    unit_cell = UnitCell(1, 1, 1, gamma=60)

    x1 = LatticeSite(0, 0, 0.5, 0, 0, 1, name="X1")
    x2 = LatticeSite(0, 0, 0, 0, 1, 0, name="X2")

    sites = [x1, x2]

    exchanges = generate_exchanges(sites=sites,
                                  unit_cell=unit_cell,
                                  max_distance=0.6,
                                  coupling_type=HeisenbergCoupling,
                                  j=-1)

    energies = []
    for supercell in [TrivialSupercell(), TrivialSupercell((2,2,2)), TrivialSupercell((1,2,3))]:

        s = Structure(sites, unit_cell=unit_cell, supercell=supercell)
        hamiltonian = Hamiltonian(s, exchanges)
        minimiser = ClassicalEnergyMinimisation(hamiltonian)

        energies.append(minimiser.energy())

    energies = np.array(energies)

    assert np.allclose(energies[1:], energies[0]), "all energies should be the same"



def test_energy_behaviour_summation_supercell():
    """ Check energy calculation is consistent for summation supercells

    Strategy, have spins pointing in different relative directions in different cells

    Check that:
    E(UU, UU) = E(DD, UU)
    E(UU, UD) = D(UU, DU) = 0
    """

    unit_cell = UnitCell(1, 1, 1, gamma=60)

    x1 = LatticeSite(0, 0, 0.5, supercell_moments=[[0, 0, 1],[0,0,0]], name="X1")
    x2 = LatticeSite(0, 0, 0, supercell_moments=[[0, 0, 0], [0, 0, 1]], name="X2")

    sites = [x1, x2]

    exchanges = [HeisenbergCoupling(x1, x2, j=-1)]

    energies = []
    energy_factors = []
    for supercell, energy in [(SummationSupercell([
                                  CommensuratePropagationVector(1,1,1),
                                  CommensuratePropagationVector(1,1,1)]), 1),
                              (SummationSupercell([
                                  CommensuratePropagationVector(1,1,1/2),
                                  CommensuratePropagationVector(1,1,1/2)]), 1),
                              (SummationSupercell([
                                  CommensuratePropagationVector(1, 1, 1/2),
                                  CommensuratePropagationVector(1, 1, 1)]), 0),
                              (SummationSupercell([
                                  CommensuratePropagationVector(1, 1, 1),
                                  CommensuratePropagationVector(1, 1, 1/2)]), 0) ]:

        s = Structure(sites, unit_cell=unit_cell, supercell=supercell)
        hamiltonian = Hamiltonian(s, exchanges)

        hamiltonian.expanded().print_summary()

        minimiser = ClassicalEnergyMinimisation(hamiltonian)

        energies.append(minimiser.energy())
        energy_factors.append(energy)

    # Scale

    ratio = energies[0] / energy_factors[0]

    expected_energies = [ratio * factor for factor in energy_factors]

    energies = np.array(energies)
    expected_energies = np.array(expected_energies)

    assert np.allclose(energies, expected_energies), "energies should match expected energies"

def test_energy_behaviour_summation_supercell():
    """ Check energy calculation is consistent for summation supercells

    Strategy, have spins pointing in different relative directions in different cells

    Check that:
    E(UU, UU) = E(DD, UU)
    E(UU, UD) = D(UU, DU) = 0
    """

    unit_cell = UnitCell(1, 1, 1, gamma=60)

    x1 = LatticeSite(0, 0, 0.5, supercell_moments=[[0, 0, 1],[0,0,0]], name="X1")
    x2 = LatticeSite(0, 0, 0, supercell_moments=[[0, 0, 0], [0, 0, 1]], name="X2")

    sites = [x1, x2]

    exchanges = [HeisenbergCoupling(x1, x2, j=-1)]

    energies = []
    energy_factors = []
    for supercell, energy in [(SummationSupercell([
                                  CommensuratePropagationVector(1,1,1),
                                  CommensuratePropagationVector(1,1,1)]), 1),
                              (SummationSupercell([
                                  CommensuratePropagationVector(1,1,1/2),
                                  CommensuratePropagationVector(1,1,1/2)]), 1),
                              (SummationSupercell([
                                  CommensuratePropagationVector(1, 1, 1/2),
                                  CommensuratePropagationVector(1, 1, 1)]), 0),
                              (SummationSupercell([
                                  CommensuratePropagationVector(1, 1, 1),
                                  CommensuratePropagationVector(1, 1, 1/2)]), 0) ]:

        s = Structure(sites, unit_cell=unit_cell, supercell=supercell)
        hamiltonian = Hamiltonian(s, exchanges)

        hamiltonian.expanded().print_summary()

        minimiser = ClassicalEnergyMinimisation(hamiltonian)

        energies.append(minimiser.energy())
        energy_factors.append(energy)

    # Scale

    ratio = energies[0] / energy_factors[0]

    expected_energies = [ratio * factor for factor in energy_factors]

    energies = np.array(energies)
    expected_energies = np.array(expected_energies)

    assert np.allclose(energies, expected_energies), "energies should match expected energies"


def test_energy_behaviour_rotation_supercell():
    """ Check energy calculation is consistent for rotation supercells

    Strategy, have spins pointing in different relative directions in different cells

    Check that:
    E(UU) = -E(UD)
    """

    unit_cell = UnitCell(1, 1, 1, gamma=60)

    x = LatticeSite(0, 0, 0.5, 0,0,1, name="X")

    exchanges = [HeisenbergCoupling(x, x, cell_offset=CellOffset(0,0,1), j=-1)]

    energies = []
    energy_factors = []
    for supercell, energy in [
        (TransformationSupercell(
            [(CommensuratePropagationVector(1,1,1), RotationTransform((1,0,0)))]), 1),
        (TransformationSupercell(
            [(CommensuratePropagationVector(1,1,1/2), RotationTransform((1,0,0)))]), -1) ]:

        s = Structure([x], unit_cell=unit_cell, supercell=supercell)
        hamiltonian = Hamiltonian(s, exchanges)

        hamiltonian.expanded().print_summary()

        minimiser = ClassicalEnergyMinimisation(hamiltonian)

        energies.append(minimiser.energy())
        energy_factors.append(energy)

    # Scale

    ratio = energies[0] / energy_factors[0]

    expected_energies = [ratio * factor for factor in energy_factors]

    energies = np.array(energies)
    expected_energies = np.array(expected_energies)

    assert np.allclose(energies, expected_energies), "energies should match expected energies"

def test_rotation_supercell_error():
    """ Test that we throw an error if called on an incommensurate structure """

    # TODO, add when integrated with Duc's changes

    # with pytest.raises(TypeError):
    #     pass

def test_optimise_transformation_supercell():
    """ Test optimisation of a supercell where moments need to be as unaligned as possible """

    a = LatticeSite(0.25,0,0, 1, 0, 0, name="A")
    x = LatticeSite(0.5, 0, 0, 1, 1, 1, name="X")

    sites = [a, x]

    pv = CommensuratePropagationVector(1,1,1/3)
    structure = Structure(sites, UnitCell(1,1,1),
                          supercell=TransformationSupercell([(pv, RotationTransform([1,0,0]))]))

    couplings = [HeisenbergCoupling(a, x, 1)]

    hamiltonian = Hamiltonian(structure, couplings)

    optimised = hamiltonian.ground_state(fixed=[a])

    x_new = optimised.sites_by_name("X")[0]


    base_moment_direction = x_new.base_moment / np.sqrt(np.sum(x_new.base_moment**2))

    assert np.isclose(np.dot(base_moment_direction, [1,0,0]), -1)

def test_optimise_summation_supercell():
    """ A system where we optimise one of the in the supercell but not the other, just check it doesn't crash """
    x = LatticeSite(0.25,0,0, supercell_moments=[[0,0,1], [0,0,0]], name="X")
    y = LatticeSite(0.75,0,0, supercell_moments=[[0,0,0], [1,1,1]], name="Y")

    couplings = [HeisenbergCoupling(x, y, j=1)]

    pv1 = CommensuratePropagationVector(1,1, 1/2)
    pv2 = CommensuratePropagationVector(1,1, 1/3, phase=np.pi/2)

    supercell = SummationSupercell([pv1, pv2])

    structure = Structure([x, y], unit_cell=UnitCell(1,1,1), supercell=supercell)

    hamiltonian = Hamiltonian(structure, couplings)

    optimised = hamiltonian.ground_state(fixed=[x])

