import numpy as np

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

