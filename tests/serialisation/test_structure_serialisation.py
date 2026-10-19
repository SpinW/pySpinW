""" Tests for structure serialisation

As we have tests for site, unit_cell, spacegroup and supercell, we only need to test that things are combined
correctly.

The main thing to check is that the symmetry implied sites are correctly serialised

"""

import pytest

from pyspinw import TiledSupercell, CommensuratePropagationVector, RotationTransform, TransformationSupercell, UnitCell, \
    LatticeSite, Structure
from pyspinw.symmetry.group import database

input_data = [(CommensuratePropagationVector(0, 0, 1 / 2), RotationTransform([0,1,0])),
                (CommensuratePropagationVector(1 / 3, 1 / 3, 1 / 3), RotationTransform([1,0,0]))]


spacegroups = database.spacegroups[::10] # Just some of them
supercells = [TiledSupercell(scaling=(2,3,4)),
              TransformationSupercell(input_data, scaling=(1, 3, 5))]
unit_cells = [UnitCell(1,1,1),
              UnitCell(2,3,3, gamma=120)]

@pytest.mark.parametrize("spacegroup", spacegroups)
@pytest.mark.parametrize("supercell", supercells)
@pytest.mark.parametrize("unit_cell", unit_cells)
def test_structure_serialisation(spacegroup, supercell, unit_cell):
    sites = [LatticeSite(1/2,0,1/2, name="S1"),
             LatticeSite(0,0,0, name="S2"),
             LatticeSite(0.1, 0.2, 0.3, name="S3")]

    # Need to skip checks on inputs, because we've been a bit careless about whether the structure is actually valid
    # Shouldn't be a problem for serialisation testing though.
    structure = Structure(sites, unit_cell=unit_cell, spacegroup=spacegroup, supercell=supercell, skip_checks=True)

    json = structure.serialise()

    deserialised = Structure.deserialise(json)

    assert isinstance(deserialised, Structure)
    assert len(structure.sites) == len(deserialised.sites)

    # Check unit cell
    assert structure.unit_cell.a == deserialised.unit_cell.a
    assert structure.unit_cell.b == deserialised.unit_cell.b
    assert structure.unit_cell.c == deserialised.unit_cell.c
    assert structure.unit_cell.alpha == deserialised.unit_cell.alpha
    assert structure.unit_cell.beta == deserialised.unit_cell.beta
    assert structure.unit_cell.gamma == deserialised.unit_cell.gamma

    # Check supercell is the same type
    assert structure.supercell.__class__ == deserialised.supercell.__class__

    # Check the spacegroup
    assert structure.spacegroup.hall_number == deserialised.spacegroup.hall_number

    # Check the summary string, this should catch whether the sites are implied or not
    assert structure.text_summary == deserialised.text_summary