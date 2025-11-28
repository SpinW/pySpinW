""" Information/structures around symmetry state """

from dataclasses import dataclass

from pyspinw.symmetry.group import MagneticSpaceGroup, SpaceGroup, magnetic_group_symbol_lookup, \
    spacegroup_symbol_lookup, magnetic_groups
from pyspinw.symmetry.unitcell import UnitCell


class SymmetrySettings:
    """ Object to hold symmetry information together """

    def __init__(self, space_group: SpaceGroup, magnetic_group: MagneticSpaceGroup, unit_cell: UnitCell):

        self._space_group = space_group
        self._magnetic_group = magnetic_group
        self._unit_cell = unit_cell

    @property
    def space_group(self):
        """ Get the space group """
        return self._space_group

    @property
    def magnetic_group(self):
        """ Get the magnetic group """
        return self._magnetic_group

    @property
    def unit_cell(self):
        """ Get the unit cell """
        return self._unit_cell

DEFAULT_SYMMETRY = SymmetrySettings(
    space_group=spacegroup_symbol_lookup["P 1"],
    magnetic_group=magnetic_group_symbol_lookup["P1.1"],
    unit_cell=UnitCell(1,1,1)
)
