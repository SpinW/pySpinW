"""Base classes for PySpinW

This is an abstract outline, the actual implementations are in different files

"""
from dataclasses import dataclass

# pylint: disable=R0903

import numpy as np

from pyspinw.cell_offsets import CellOffset, CellOffsetCoercible
from pyspinw.checks import check_sizes
from pyspinw.serialisation import SPWSerialisable, SPWSerialisationContext, SPWDeserialisationContexGroup, \
    SPWDeserialisationContext
from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell


class MagneticStructure(SPWSerialisable):
    """Base class for representations of the Magnetic Structures"""

    serialisation_name = "structure"

    def __init__(self):
        pass

    def _serialise(self, context: SPWSerialisationContext):
        return {}

    @staticmethod
    def _deserialise(data: dict, context: SPWDeserialisationContext):
        return MagneticStructure()



class Data:
    """Placeholder"""

    def __init__(self, data):
        self.data = data

    @property
    def q(self) -> np.ndarray:
        raise NotImplementedError
