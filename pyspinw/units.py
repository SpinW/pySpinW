""" Units """
from enum import Enum

class CoordsUnits(Enum):
    """ Types of coordinate system """

    XYZ = 'xyz' # Cartesian
    LU = 'lu'   # Lattice units


class IntensityUnits(Enum):
    """ Output intensity units """

    PERCELL = 'spin length per cell'
    PERSPIN = 'spin length per spin'
    BARN = 'barns / sr / meV / atom'


def intensity_units(name):
    """ Creates an intensity unit from a string name """
    if not isinstance(name, str):
        return IntensityUnits(name)
    elif 'cell' in name:
        return IntensityUnits.PERCELL
    elif 'spin' in name:
        return IntensityUnits.PERSPIN
    elif 'barn' in name:
        return IntensityUnits.BARN
