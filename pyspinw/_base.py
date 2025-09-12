"""Base classes for PySpinW

This is an abstract outline, the actual implementations are in different files

"""

# pylint: disable=R0903

import numpy as np

from pyspinw.cell_offsets import CellOffset, CellOffsetCoercible
from pyspinw.checks import check_sizes
from pyspinw.serialisation import SPWSerialisable, SPWSerialisationContext, SPWDeserialisationContexGroup, \
    SPWDeserialisationContext
from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell


class Data:
    """Placeholder"""

    def __init__(self, data):
        self.data = data

    @property
    def q(self) -> np.ndarray:
        raise NotImplementedError
