import pytest

from pyspinw.symmetry.group import database, SpaceGroup
from pyspinw.symmetry.system import lattice_systems, Cubic, Hexagonal, Rhombohedral, Tetragonal, Orthorhombic, \
    Monoclinic, Triclinic, IncompatibleUnitCell
from pyspinw.symmetry.unitcell import UnitCell


def test_all_systems_present():
    """ Test that all the Bravais lattices are there are correct """
    expected_names = sorted([
        "aP",
        "mP", "mS",
        "oP", "oS", "oI", "oF",
        "tP", "tI",
        "hR",
        "hP",
        "cP", "cI", "cF"
    ])

    assert len(expected_names) == 14 # Sanity check

    names = []
    for lattice_system in lattice_systems:
        names += [lattice_system.letter + bravais.letter for bravais in lattice_system.bravais_options.bravias]

    names.sort()

    for expected, calculated in zip(expected_names, names):
        assert expected == calculated

#
# Tests for unit cell validation
#

examples = [
    (Triclinic(), UnitCell(1,2,3, 60, 70, 80)),
    (Monoclinic(), UnitCell(1,2,3, 90, 90, 80)),
    (Orthorhombic(), UnitCell(1,2,3, 90,90,90)),
    (Tetragonal(), UnitCell(1,1,3, 90,90,90)),
    (Rhombohedral(), UnitCell(7,7,7,50, 50, 50)),
    (Hexagonal(), UnitCell(1,1,3, 90, 90, 120)),
    (Cubic(), UnitCell(7,7,7,90,90,90))]

@pytest.mark.parametrize("system_and_cell", examples)
@pytest.mark.parametrize("other_system", [system for system, _ in examples])
def test_unit_cell_validation(system_and_cell, other_system):
    """ LatticeSystem.validate should raise an error iff the unit cell parameters don't fit"""
    system, cell = system_and_cell
    if system.name == other_system.name:
        # If it is the right system validate should not throw an error
        other_system.validate(cell)
    else:
        with pytest.raises(IncompatibleUnitCell):
            # If it is the wrong system it should throw an error
            other_system.validate(cell)


# Build a list of spacegroups with different lattice systems
example_groups: list[SpaceGroup] = []
for lattice_system in lattice_systems:
    name = lattice_system.name

    for spacegroup in database.spacegroups:
        if spacegroup.lattice_system.name == name:
            example_groups.append(spacegroup)
            break


example_unit_cell_access: dict[str, dict[str, float]] = {
    Triclinic.name: {"a": 1, "b": 2, "c": 3, "alpha": 40, "beta": 50, "gamma": 60},
    Monoclinic.name: {"a": 1, "b": 2, "c": 3, "gamma": 60},
    Orthorhombic.name: {"a": 1, "b": 2, "c": 3},
    Tetragonal.name: {"a": 1, "c": 2},
    Rhombohedral.name: {"a": 7, "alpha": 60},
    Hexagonal.name: {"a": 1, "c": 3},
    Cubic.name: {"a": 7}
}

@pytest.mark.parametrize("group", example_groups)
def test_create_unit_cell(group: SpaceGroup):
    """ Check that providing the right parameters to different space group's `create_unit_cell` works """
    params = example_unit_cell_access[group.lattice_system.name]

    group.create_unit_cell(**params)
