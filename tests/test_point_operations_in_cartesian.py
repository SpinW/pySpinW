import numpy as np
import pytest

from pyspinw import UnitCell
from pyspinw.symmetry.group import database, SpaceGroup

def labeller(sg: SpaceGroup):
    """ Label for spacegroups to go on the test names"""
    if sg.choice is None:
        return f"{sg.symbol}, ({sg.hall_number}, {sg.lattice_system.name})"
    else:
        return f"{sg.symbol} [{sg.choice}], ({sg.hall_number}, {sg.lattice_system.name})"


@pytest.mark.parametrize("sg", database.spacegroups, ids=labeller)
@pytest.mark.parametrize("cell", [
    UnitCell(1,1,1),
    UnitCell(3,4,5, 57, 59, 110)],ids=["Symmetric Cell", "Asymmetric Cell"])
def test_action_of_spacegroup_on_spins(sg: SpaceGroup, cell: UnitCell):
    """ Checks that the spacegroup acting on spins is"""
    test_cell = sg.lattice_system.constrain(cell)

    for operation in sg:
        op = operation.point_operation_for_spins(test_cell)

        assert np.allclose(op @ op.T, np.eye(3)), f"Failed on {sg}: {operation}"

