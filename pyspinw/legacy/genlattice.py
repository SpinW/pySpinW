""" This will hold code that implements genlattice for pyspinw"""
from dataclasses import dataclass

import numpy as np
from numpy._typing import ArrayLike

from pyspinw.checks import check_sizes
from pyspinw.symmetry.group import SpaceGroup
from pyspinw.interface import spacegroup as lookup_spacegroup
from pyspinw.symmetry.unitcell import UnitCell

@dataclass
class GenLatticeResult:
    """ Result from a genlattice call """

    spacegroup: SpaceGroup
    unit_cell: UnitCell



@check_sizes(angle=(3,), angled=(3,), lat_const=(3,), allow_nones=True, force_numpy=True)
def genlattice(
        lat_const: ArrayLike,
        angle: ArrayLike | None = None,
        angled: ArrayLike | None = None,
        spgr: SpaceGroup | str = "P1"):
    """ pySpinW version of MATLAB spinW's genlattice

    This works similarly, though the definition of spacegroups is a little different.

    :param lat_const: lattice constants a,b and c
    :param angle: unit cell angles alpha, beta and gamma in radians (you need to specify either this or angled)
    :param angled: unit cell angles alpha, beta and gamma in degrees (you need to specify either this or angle)
    :param spgr: spacegroup, either a search string (see `pyspinw.interface.spacegroup`) or a SpaceGroup object

    :return: GenLatticeResult, an object that holds the information generated, which can be fed to genmagstr
    """
    if isinstance(spgr, str):
        spacegroup = lookup_spacegroup(spgr)

    elif isinstance(spgr, SpaceGroup):
        spacegroup = spgr

    else:
        raise TypeError("Expected `spgr` type to be either `str` or `SpaceGroup`")


    if angle is None and angled is None:
        raise ValueError("Expected to have either 'angle' or 'angled' specified.")

    if angle is not None and angled is not None:
        raise ValueError("Expected either 'angle' or 'angled' to be specified, not both.")

    if angled is None:
        angled = (180/np.pi)*angle

    # Check compatability with cell definition

    base_cell = UnitCell(a=lat_const[0],
                         b=lat_const[1],
                         c=lat_const[2],
                         alpha=angled[0],
                         beta=angled[1],
                         gamma=angled[2])

    # Validate the cell
    spacegroup.lattice_system.validate(base_cell)

    return GenLatticeResult(spacegroup, base_cell)

if __name__ == "__main__":
    res = genlattice(angled=(40,50,60), lat_const=(1,2,3), spgr="p1")
    print(res)
