""" Information/structures around symmetry state """

from dataclasses import dataclass

from pyspinw.symmetry.group import MagneticSpaceGroup, SpaceGroup, magnetic_group_symbol_lookup, \
    spacegroup_symbol_lookup, magnetic_groups
from pyspinw.symmetry.unitcell import UnitCell


@dataclass
class SymmetrySettings:
    """ Object to hold symmetry information together """

    space_group: SpaceGroup
    magnetic_group: MagneticSpaceGroup
    unit_cell: UnitCell

DEFAULT_SYMMETRY = SymmetrySettings(
    space_group=spacegroup_symbol_lookup["P 1"],
    magnetic_group=magnetic_group_symbol_lookup["P1.1"],
    unit_cell=UnitCell(1,1,1)
)
