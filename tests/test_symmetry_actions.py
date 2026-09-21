import numpy as np


import pytest

from pyspinw import HeisenbergExchange
from symmetry.exchange_test_cases import heisenberg_cases, diagonal_cases

q_points = np.concatenate((
        np.linspace(-0.9, -0.1, 9),
        np.linspace(0.1, 0.9, 9)))

qx, qy, qz = np.meshgrid(q_points, q_points, q_points)

qs = np.concatenate((
    qx.reshape(-1, 1),
    qy.reshape(-1, 1),
    qz.reshape(-1, 1)), axis=1)


@pytest.mark.parametrize("hamiltonian", heisenberg_cases + diagonal_cases)
def test_tests_for_symmetry_filled(hamiltonian):
    """ This makes sure the tests we have should fail if things are wrong.
    Make an asymmetric system, then apply operations to the system *but not to q*,
    it should give different spinwave results"""

    initial_energies, _ = hamiltonian._energies_and_intensities(qs, use_rotating=False)

    all_not_failing = True

    for op in hamiltonian.structure.spacegroup:
        transformed_hamiltonian = hamiltonian.symmetry_transformed(op)

        compare_energies, _ = transformed_hamiltonian._energies_and_intensities(qs, use_rotating=False)

        transformed_hamiltonian.print_summary()
        for exchange in transformed_hamiltonian.exchanges:
            print(exchange.exchange_matrix)

        print(np.max(np.abs(np.array(compare_energies) - initial_energies)))

        all_not_failing = all_not_failing and np.allclose(initial_energies, compare_energies)

    assert not all_not_failing, "A good test Hamiltonian should not be symmetry invariant"


@pytest.mark.skip(reason="It's unclear whether this test is correctly specified.")
@pytest.mark.parametrize("hamiltonian", heisenberg_cases)
def test_hamiltonian_symmetry_operations_apply_correctly(hamiltonian):

    """ Make an asymmetric system, then apply operations to the system and q, it should give the same spinwave result"""


    # qs_in_lattice = hamiltonian.structure.unit_cell.spin_cartesian_to_lattice_units(qs)
    qs_in_lattice = hamiltonian.structure.unit_cell.spin_lattice_units_to_cartesian(qs)
    initial_energies, _ = hamiltonian._energies_and_intensities(qs_in_lattice, use_rotating=False)


    for op in hamiltonian.structure.spacegroup:
        transform = op.point_operation_in_cartesian(hamiltonian.structure.unit_cell)
        transformed_qs = (transform @ qs.T).T

        # transformed_qs_in_lattice = hamiltonian.structure.unit_cell.spin_cartesian_to_lattice_units(transformed_qs)
        transformed_qs_in_lattice = hamiltonian.structure.unit_cell.spin_lattice_units_to_cartesian(transformed_qs)


        transformed_hamiltonian = hamiltonian.symmetry_transformed(op)

        compare_energies, _ = transformed_hamiltonian._energies_and_intensities(
            transformed_qs_in_lattice, use_rotating=False)


        # import matplotlib.pyplot as plt
        # plt.plot(np.arange(len(qs)), initial_energies, color='k')
        # plt.plot(np.arange(len(qs)), compare_energies, color='r')
        # plt.show()


        max_difference = np.max(np.abs(np.array(initial_energies) - compare_energies))

        close = max_difference < 1e-12

        assert close, f"Failed on operation {op}, max difference {max_difference}"

def sort_energies(energies: list[np.ndarray]):
    """ Sort energies """
    return [np.sort(energy) for energy in energies]

@pytest.mark.parametrize("hamiltonian", heisenberg_cases + diagonal_cases)
def test_symmetry_filled_symmetric_exchange_cases(hamiltonian):

    """ Make an asymmetric system, fill it by symmetry, then apply operations to the system *but not to q*,
    it should give the same spinwave result"""



    hamiltonian = hamiltonian.symmetry_filled() # Get hamiltonian that should obey symmetry

    initial_energies, _ = hamiltonian._energies_and_intensities(qs, use_rotating=False)
    initial_energies = sort_energies(initial_energies)

    for op in hamiltonian.structure.spacegroup:
        transformed_hamiltonian = hamiltonian.symmetry_transformed(op)

        compare_energies, _ = transformed_hamiltonian._energies_and_intensities(qs, use_rotating=False)
        compare_energies = sort_energies(compare_energies)

        close = np.allclose(initial_energies, compare_energies)

        # import matplotlib.pyplot as plt
        # plt.plot(np.arange(len(qs)), initial_energies, color='k')
        # plt.plot(np.arange(len(qs)), compare_energies, color='r')
        # plt.show()

        assert close, f"Failed on operation {op}"


@pytest.mark.parametrize("hamiltonian", heisenberg_cases)
def test_fill_symmetry_makes_heisenbergs_from_heisenbergs(hamiltonian):

    """ Do symmetry fill on an asymmetric system with Heisenberg exchanges, the result should be Heisenbergs too """

    filled = hamiltonian.symmetry_filled()

    for exchange in filled.exchanges:
        assert isinstance(exchange, HeisenbergExchange), ("Heisenberg exchanges should transform to Heisenberg "
                                                          f"exchanges, but got {exchange} (exchange matrix is"
                                                          f"{exchange.exchange_matrix})")



