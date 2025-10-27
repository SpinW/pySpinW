import pytest

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
