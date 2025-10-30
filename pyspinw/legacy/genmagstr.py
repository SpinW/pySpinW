""" Implementation of python version of genmagstr """
from dataclasses import dataclass
from enum import Enum

from numpy._typing import ArrayLike

from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell


class GenMagStrMode(Enum):
    """ Different ways of generating magnetic structures using genmagstr"""

    RANDOM = 'random'
    DIRECT = 'direct'
    TILE = 'tile'
    HELICAL = 'helical'
    ROTATE = 'rotate'
    FUNC = 'func' # Let's not bother with this though
    FOURIER = 'fourier'
    EXTEND = 'extend' # Deprecated

class UnitSystem(Enum):
    """ Types of coordinate system """

    XYZ = 'xyz' # Cartesian
    LU = 'lu'   # Lattice units

@dataclass
class GenMagStrResult:
    sites: list[LatticeSite]


def genmagstr(
        mode: GenMagStrMode | str,
        sites: ArrayLike | list[LatticeSite],
        s: ArrayLike,
        phi: float = 0.0,
        phi_d: float = 0.0,
        n_ext: tuple[int, int, int] = (1,1,1),
        k: tuple[float, float, float] = (0,0,0) | None,
        n: tuple[float, float, float] = (0,0,1),
        unit: UnitSystem | str = UnitSystem.XYZ,
        unit_cell: UnitCell = UnitCell(1,1,1),
        norm: bool = True,
        r0: bool = True):
    """ TODO: Needs a docstring"""
    # TODO: Make this apply to an object that contains unit cell and sites,
    #  will depend on the output from other functions used to create sites (i.e. genlattice)

    # Convert moments out of unit system to xyz
    match unit:
        case UnitSystem.XYZ:
            moments = s
        case UnitSystem.LU:
            moments = unit_cell.fractional_to_cartesian(s)


    match mode:
        case GenMagStrMode.EXTEND:
            raise NotImplementedError("'extend' mode has been deprecated, use 'tile' instead")

        case GenMagStrMode.FUNC:
            raise NotImplementedError("'func' mode is not implemented in pySpinW, for better control over"
                                      " structures use the python interface rather than this legacy interface")

        case GenMagStrMode.TILE:
            # This is just the idenity supercell
            pass

        case GenMagStrMode.DIRECT:
            pass

        case GenMagStrMode.ROTATE:
            pass
