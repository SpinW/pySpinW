""" This will hold code that implements genlattice for pyspinw"""
from pyspinw.checks import check_sizes
from pyspinw.site import LatticeSite
from pyspinw.symmetry.group import SpaceGroup
from pyspinw.symmetry.unitcell import UnitCell


class GenLatticeResult:
    spacegroup: SpaceGroup
    unit_cell: UnitCell

@check_sizes()
def genlattice(angles: ):
    pass