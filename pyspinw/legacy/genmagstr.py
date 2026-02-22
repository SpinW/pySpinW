""" Implementation of python version of genmagstr """
from enum import Enum
import numpy as np
from numpy._typing import ArrayLike
from typing import Callable

from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structures import Structure
from pyspinw.symmetry.supercell import TrivialSupercell, SummationSupercell, CommensuratePropagationVector


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


def genmagstr(
        sites: ArrayLike | list[LatticeSite],
        unit_cell: UnitCell,
        mode: GenMagStrMode | str,
        phi: ArrayLike | None = None,
        phid: ArrayLike | None = None,
        nExt: ArrayLike | None = None,
        k: ArrayLike | None = None,
        n: ArrayLike | None = None,
        S: ArrayLike | None = None,
        unit: UnitSystem | str = UnitSystem.XYZ,
        epsilon: float = 1e-5,
        func: Callable | None = None,
        x0: ArrayLike | None = None,
        norm: bool = True,
        r0: bool = True) -> Structure:
    """ generates a magnetic structure using Matlab-like syntax """
    mode, unit = (GenMagStrMode(mode), UnitSystem(unit))

    # Some parameter checking
    nExt = (1,1,1) if nExt is None else tuple(*nExt)

    # Convert moments out of unit system to xyz
    if S is not None and unit == UnitSystem.LU:
        S = unit_cell.fractional_to_cartesian(S)

    match mode:
        case GenMagStrMode.TILE | GenMagStrMode.EXTEND:
            if unit == UnitSystem.LU:
                norm_transform = unit_cell._xyz / np.sqrt(np.sum(unit_cell._xyz**2, axis=1)).reshape(-1, 1)
                for site in sites:
                    Stmp = site._base_moment @ norm_transform
                    site._base_moment = Stmp * site._magnitude / np.linalg.norm(Stmp)
                    site._moment_data = np.array([site._base_moment])
            return Structure(sites, unit_cell, supercell=TrivialSupercell(scaling=nExt))

        case GenMagStrMode.RANDOM:
            pass

        case GenMagStrMode.HELICAL:
            # Note: not all cases handled by Matlab are handled by this
            if n is None or k is None:
                raise RuntimeError('You must provide the rotation axis "n" and propagation vector "k"')
            n = np.array(n, dtype='double')
            norm_transform = unit_cell._xyz / np.sqrt(np.sum(unit_cell._xyz**2, axis=1)).reshape(-1, 1)
            # Make the basis complex
            for ii, site in enumerate(sites):
                if S is not None:
                    S = S * site._magnitude / np.linalg.norm(S)
                    site._base_moment = S[:,ii] + 1j * np.cross(n, S[:,ii])
                elif unit == UnitSystem.LU:
                    Stmp = site._base_moment @ norm_transform
                    Stmp = Stmp * site._magnitude / np.linalg.norm(Stmp)
                    site._base_moment = Stmp + 1j * np.cross(n, Stmp)
                else:
                    site._base_moment = site._base_moment * site._magnitude / np.linalg.norm(site._base_moment)
                    site._base_moment = site._base_moment + 1j * np.cross(n, site._base_moment)
                site._moment_data = np.array([site._base_moment])
            k = CommensuratePropagationVector(k[0], k[1], k[2])
            return Structure(sites, unit_cell, supercell=SummationSupercell(propagation_vectors=[k]))

        case GenMagStrMode.DIRECT:
            raise NotImplementedError("'direct' mode to be implemented")

        case GenMagStrMode.ROTATE:
            raise NotImplementedError("'rotate' mode to be implemented")

        case GenMagStrMode.FUNC:
            raise NotImplementedError("'func' mode is not implemented in pySpinW, for better control over"
                                      " structures use the python interface rather than this legacy interface")

