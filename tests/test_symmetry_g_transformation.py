import numpy as np
import pytest
from numpy._typing import ArrayLike

from pyspinw import UnitCell, LatticeSite
from pyspinw.symmetry.group import database, SpaceGroup

def labeller(sg: SpaceGroup):
    """ Label for spacegroups to go on the test names"""
    if sg.choice is None:
        return f"{sg.symbol}, ({sg.hall_number}, {sg.lattice_system.name})"
    else:
        return f"{sg.symbol} [{sg.choice}], ({sg.hall_number}, {sg.lattice_system.name})"

magnetic_vectors = [
    [2,0,0],
    [0,1,1],
    [-3,-2,-7]
]

spins = [[1,0,0],[0,1,0],[0,0,1],[1,1,0],[0,1,1],[1,1,1]]

m = np.array([[1,2,3],[4,5,6],[7,8,9]])
g_tensors = [
    np.eye(3),
    2*np.eye(3),
    np.diag([2,3,4]),
    m + m.T]

sites = [
    LatticeSite(0,0,0,supercell_spins=spin, g=g_tensor)
    for spin in spins
    for g_tensor in g_tensors
]

@pytest.mark.parametrize("sg", database.spacegroups, ids=labeller)
@pytest.mark.parametrize("cell", [
    UnitCell(1,1,1),
    UnitCell(3,4,5, 57, 59, 110)],ids=["Most Symmetric", "Least Symmetric"])
@pytest.mark.parametrize("magnetic_vector", magnetic_vectors)
def test_space_operations_on_spins_preserves_length(
        sg: SpaceGroup,
        cell: UnitCell,
        magnetic_vector: ArrayLike):
    """ Checks that the spacegroup acting on spins preserves length"""
    test_cell = sg.lattice_system.constrain(cell)

    for site in sites:
        original_zeeman = magnetic_vector @ site.g @ site.spin

        for operation in sg:
            new_site = site.symmetry_transformed(operation, test_cell)

            transform = operation.point_operation_in_cartesian(test_cell)
            new_magnetic_vector = transform @ magnetic_vector

            zeeman = new_magnetic_vector @ new_site.g @ new_site.spin

            assert np.isclose(zeeman, original_zeeman), f"Failed on {site}, {sg}: {operation}"

