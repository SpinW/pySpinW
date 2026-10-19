""" Check when building a structure warns that the unit cell is more regular than its lattice system requires """

import logging

from pyspinw.site import LatticeSite
from pyspinw.structure import Structure
from pyspinw.symmetry.group import database
from pyspinw.symmetry.unitcell import UnitCell


def _site():
    return LatticeSite(0, 0, 0, 0, 0, 1, name="S")


def test_default_spacegroup_does_not_warn_about_the_cell(caplog):
    """ Any cell is a valid triclinic cell, so P1 (the default) has nothing to warn about """
    with caplog.at_level(logging.WARNING, logger="Structure"):
        structure = Structure([_site()], UnitCell(3.0, 8.0, 8.0))
        structure.without_nonmagnetic()

    assert caplog.records == []


def test_orthorhombic_spacegroup_with_a_cubic_cell_still_warns(caplog):
    """ A cubic cell is more regular than an orthorhombic space group requires """
    with caplog.at_level(logging.WARNING, logger="Structure"):
        Structure([_site()], UnitCell(5.0, 5.0, 5.0), spacegroup=database.spacegroup_by_name("Pmmm"))

    assert len(caplog.records) == 1
    assert "Orthorhombic" in caplog.records[0].getMessage()
