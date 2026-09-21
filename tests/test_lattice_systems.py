import pytest

from pyspinw.symmetry.group import database, SpaceGroup
from pyspinw.symmetry.system import lattice_systems, Cubic, Hexagonal, Rhombohedral, Tetragonal, Orthorhombic, \
    MonoclinicA, MonoclinicB, MonoclinicC, Triclinic, IncompatibleUnitCell
from pyspinw.symmetry.unitcell import UnitCell


lattice_systems = [lattice for lattice in lattice_systems if len(lattice.letter) < 2]

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
