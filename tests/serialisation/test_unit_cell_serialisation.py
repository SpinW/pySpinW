""" Tests for serialisation of UnitCell """

import pytest

import numpy as np

from pyspinw import UnitCell

cells = [
    UnitCell(1,1,1),
    UnitCell(2,3,4),
    UnitCell(1.234523452345,2.234523452,3.234523452345,
             87.2345234523425,88.2345234523452345,89.324523453245),
    UnitCell(1,1,1,ab_normal=(1,2,3), direction=(4,5,6))
]

@pytest.mark.parametrize("cell", cells)
def test_unit_cell_serialisation(cell):
    """ Check that unit cells serialise correctly """

    json = cell.serialise()
    deserialised = UnitCell.deserialise(json)

    assert isinstance(deserialised, UnitCell)

    assert cell.a == deserialised.a
    assert cell.b == deserialised.b
    assert cell.c == deserialised.c
    assert cell.alpha == deserialised.alpha
    assert cell.beta == deserialised.beta
    assert cell.gamma == deserialised.gamma

    assert cell.direction == deserialised.direction
    assert cell.ab_normal == deserialised.ab_normal

    assert np.all(cell._xyz == deserialised._xyz)
    assert np.all(cell._xyz_spins == deserialised._xyz_spins)
    assert np.all(cell._xyz_inv == deserialised._xyz_inv)
    assert np.all(cell._xyz_spins_inv == deserialised._xyz_spins_inv)
