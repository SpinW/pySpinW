import numpy as np

from pyspinw import *
import pytest


test_spacegroups = [spacegroup("pmmm"), spacegroup("p4")]
test_unit_cells = [UnitCell(1,1,1), UnitCell(2,3,4), UnitCell(2,3, 4, alpha=50, beta=60)]

q_points = np.concatenate((
        np.linspace(-0.8, -0.2, 4),
        np.linspace(0.2, 0.8, 4)))

qx, qy, qz = np.meshgrid(q_points, q_points, q_points)

qs = np.concatenate((
    qx.reshape(-1, 1),
    qy.reshape(-1, 1),
    qz.reshape(-1, 1)), axis=1)

@pytest.mark.parametrize("sg", test_spacegroups)
@pytest.mark.parametrize("cell", test_unit_cells)
def test_hamiltonian_symmetry_operations_apply_correctly_to_asymmetric_exchange_cases(sg, cell):

    """ Make an asymmetric system, then apply operations to the system and q, it should give the same spinwave result"""

    initial_site = LatticeSite(0.2, 0.1, 0.3, sx=0, sy=0, sz=1/2, name="S")
    structure = Structure([initial_site], unit_cell=cell, spacegroup=sg)
    other_site, offset = structure.one_neighbour(initial_site)
    exchange = HeisenbergExchange(initial_site, other_site, j=-1.0, cell_offset=offset)
    hamiltonian = Hamiltonian(structure, [exchange])

    initial_energies, _ = hamiltonian._energies_and_intensities(qs, use_rotating=False)

    for op in sg:
        transformed_qs = (op.point_operation_matrix @ qs.T).T
        transformed_hamiltonian = hamiltonian.symmetry_transformed(op)

        compare_energies, _ = transformed_hamiltonian._energies_and_intensities(transformed_qs, use_rotating=False)

        assert np.allclose(initial_energies, compare_energies)

