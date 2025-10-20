""" This will hold code that implements genlattice for pyspinw"""
import numpy as np

from pyspinw.checks import check_sizes
from pyspinw.site import LatticeSite
from pyspinw.symmetry.group import SpaceGroup
from pyspinw.symmetry.unitcell import UnitCell


class GenLatticeResult:
    spacegroup: SpaceGroup
    unit_cell: UnitCell

@check_sizes(angle=(3,), angled=(3,), )
def genlattice(angle: np.ndarray | None = None, angled: np.ndarray | None = None, ):
    pass


