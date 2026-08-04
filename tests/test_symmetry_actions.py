import numpy as np

from pyspinw import *
import pytest

from symmetry.symmetry_test_cases import case_1, case_2

cases = [case_1(), case_2()]

q_points = np.concatenate((
        np.linspace(-0.8, -0.2, 4),
        np.linspace(0.2, 0.8, 4)))

qx, qy, qz = np.meshgrid(q_points, q_points, q_points)

qs = np.concatenate((
    qx.reshape(-1, 1),
    qy.reshape(-1, 1),
    qz.reshape(-1, 1)), axis=1)


@pytest.mark.parametrize("hamiltonian", cases)
def test_tests_for_symmetry_filled_symmetric_exchange_cases(hamiltonian):
    """ This makes sure the tests we have should fail if things are wrong.
    Make an asymmetric system, fill it, then apply operations to the system *but not to q*,
    it should give the same spinwave result"""

    initial_energies, _ = hamiltonian._energies_and_intensities(qs, use_rotating=False)

    all_not_failing = True

    for op in hamiltonian.structure.spacegroup:
        transformed_hamiltonian = hamiltonian.symmetry_transformed(op)

        compare_energies, _ = transformed_hamiltonian._energies_and_intensities(qs, use_rotating=False)

        all_not_failing = all_not_failing and np.allclose(initial_energies, compare_energies)

    assert not all_not_failing, "A good test Hamiltonian should not be symmetry invariant"



@pytest.mark.parametrize("hamiltonian", cases)
def test_hamiltonian_symmetry_operations_apply_correctly_to_symmetric_exchange_cases(hamiltonian):

    """ Make an asymmetric system, then apply operations to the system and q, it should give the same spinwave result"""


    initial_energies, _ = hamiltonian._energies_and_intensities(qs, use_rotating=False)

    for op in hamiltonian.structure.spacegroup:
        transformed_qs = (op.point_operation_matrix @ qs.T).T
        transformed_hamiltonian = hamiltonian.symmetry_transformed(op)

        compare_energies, _ = transformed_hamiltonian._energies_and_intensities(transformed_qs, use_rotating=False)

        print(transformed_hamiltonian, )

        close = np.allclose(initial_energies, compare_energies)

        assert close, f"Failed on operation {op}"

def test_g_transform():
    """ Check that the g-tensor transforms correctly"""

    # TODO

@pytest.mark.parametrize("hamiltonian", cases)
def test_symmetry_filled_symmetric_exchange_cases(hamiltonian):

    """ Make an asymmetric system, fill it, then apply operations to the system *but not to q*,
    it should give the same spinwave result"""

    hamiltonian = hamiltonian.symmetry_filled() # Get hamiltonian that should obey symmetry

    initial_energies, _ = hamiltonian._energies_and_intensities(qs, use_rotating=False)

    for op in hamiltonian.structure.spacegroup:
        transformed_hamiltonian = hamiltonian.symmetry_transformed(op)

        compare_energies, _ = transformed_hamiltonian._energies_and_intensities(qs, use_rotating=False)

        close = np.allclose(initial_energies, compare_energies)

        assert close, f"Failed on operation {op}"

