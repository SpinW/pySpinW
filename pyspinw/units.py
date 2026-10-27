"""Units."""
from enum import Enum

class CoordsUnits(Enum):
    """Types of coordinate system."""

    XYZ = 'xyz' # Cartesian
    LU = 'lu'   # Lattice units


class IntensityUnits(Enum):
    """Output intensity units."""

    PERCELL = 'spin length per cell'
    PERSPIN = 'spin length per spin'
    BARNPERATOM = 'barns / sr / meV / atom'
    BARNPERCELL = 'barns / sr / meV / cell'


def intensity_units(name):
    """Create an intensity unit from a string name.

    Parameters
    ----------
    name : IntensityUnits or str
        Intensity unit enum member or string specification.

    Returns
    -------
    IntensityUnits or None
        Resolved ``IntensityUnits`` enum member.
    """
    if not isinstance(name, str):
        return IntensityUnits(name)
    elif 'cell' in name:
        return IntensityUnits.PERCELL
    elif 'spin' in name:
        return IntensityUnits.PERSPIN
    elif 'barn' in name:
        if 'cell' in name:
            return IntensityUnits.BARNPERCELL
        else:
            return IntensityUnits.BARNPERATOM
