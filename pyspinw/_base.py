"""Base classes for PySpinW

This is an abstract outline, the actual implementations are in different files

"""
from dataclasses import dataclass

# pylint: disable=R0903

import numpy as np

from pyspinw.cell_offsets import CellOffset, CellOffsetCoercible
from pyspinw.checks import check_sizes
from pyspinw.serialisation import SPWSerialisable, SPWSerialisationContext
from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell


class MagneticStructure(SPWSerialisable):
    """Base class for representations of the Magnetic Structures"""

    def __init__(self):
        pass

    def serialise(self):
        return {}

    def deserialise(data: dict):
        return MagneticStructure()



class Data:
    """Placeholder"""

    def __init__(self, data):
        self.data = data

    @property
    def q(self) -> np.ndarray:
        raise NotImplementedError
