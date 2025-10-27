""" This will hold code that implements genlattice for pyspinw"""
import numpy as np
from numpy._typing import ArrayLike

from pyspinw.checks import check_sizes
from pyspinw.site import LatticeSite
from pyspinw.symmetry.group import SpaceGroup, spacegroup_by_name
from pyspinw.symmetry.unitcell import UnitCell
from ase.geometry.cell import cell_to_cellpar


class GenLatticeResult:
    spacegroup: SpaceGroup
    unit_cell: UnitCell



@check_sizes(angle=(3,), angled=(3,), lat_const=(3,), allow_nones=True, force_numpy=True)
def genlattice(
        lat_const: ArrayLike,
        angle: ArrayLike | None = None,
        angled: ArrayLike | None = None,
        spgr: SpaceGroup | str = "P1"):

    # a better name
    spacegroup = spgr

    if angle is None and angled is None:
        raise ValueError("Expected to have either 'angle' or 'angled' specified.")

    if angle is not None and angled is not None:
        raise ValueError("Expected either 'angle' or 'angled' to be specified, not both.")

    if angled is None:
        angled = (180/np.pi)*angle

    # Resolve the spacegroup if specified by name

    if isinstance(spacegroup, str):
        spacegroup = spacegroup_by_name(spacegroup)

    # Check compatability with cell definition

    base_cell = UnitCell(a=lat_const[0],
                         b=lat_const[1],
                         c=lat_const[2],
                         alpha=angled[0],
                         beta=angled[1],
                         gamma=angled[2])

    spacegroup.lattice_system.constrain()

    return spacegroup

if __name__ == "__main__":
    genlattice(angle=(1,2,3), lat_const=(1,2,3), spgr="p7")