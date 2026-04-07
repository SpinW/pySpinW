""" Implementation of python version of genmagstr """
from enum import Enum
import numpy as np
from numpy._typing import ArrayLike
from typing import Callable

from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structures import Structure
from pyspinw.symmetry.supercell import (TiledSupercell, RotationSupercell, SummationSupercell,
                                        CommensuratePropagationVector)


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
        magnitude: ArrayLike | None = None,
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
    if S is not None:
        if magnitude is None:
            magnitude = [np.linalg.norm(S[i,:]) for i in S.shape[0]]
        if unit == UnitSystem.LU:
            S = unit_cell.lattice_units_to_cartesian(S)
    elif magnitude is None:
        magnitude = [np.linalg.norm(site._base_moment) for site in sites]

    def _convert_moments():
        norm_transform = unit_cell._xyz / np.sqrt(np.sum(unit_cell._xyz**2, axis=1)).reshape(-1, 1)
        for i, site in enumerate(sites):
            if S is not None:
                site._base_moment = S[i,:] * magnitude[i] / np.linalg.norm(S[i,:])
            elif unit == UnitSystem.LU:
                Stmp = site._base_moment @ norm_transform
                site._base_moment = Stmp * magnitude[i] / np.linalg.norm(Stmp)
            site._moment_data = np.array([site._base_moment])

    match mode:
        case GenMagStrMode.TILE | GenMagStrMode.EXTEND:
            if S is not None or unit == UnitSystem.LU:
                _convert_moments()
            return Structure(sites, unit_cell, supercell=TiledSupercell(scaling=nExt))

        case GenMagStrMode.RANDOM:
            pass

        case GenMagStrMode.HELICAL:
            # Note: not all cases handled by Matlab are handled by this
            if n is None or k is None:
                raise RuntimeError('You must provide the rotation axis "n" and propagation vector "k"')
            if S is not None or unit == UnitSystem.LU:
                _convert_moments()
            n = np.array(n, dtype='double')
            return Structure(sites, unit_cell, supercell=RotationSupercell(perpendicular=n, propagation_vector=k))

        case GenMagStrMode.DIRECT:
            raise NotImplementedError("'direct' mode to be implemented")

        case GenMagStrMode.ROTATE:
            raise NotImplementedError("'rotate' mode to be implemented")

        case GenMagStrMode.FUNC:
            raise NotImplementedError("'func' mode is not implemented in pySpinW, for better control over"
                                      " structures use the python interface rather than this legacy interface")

