from pyspinw.gui.symmetry_settings import SymmetrySettings
from pyspinw.serialisation import SPWSerialisable, SPWSerialisationContext, SPWDeserialisationContext
from pyspinw.site import LatticeSite
from pyspinw.supercell import Supercell
from pyspinw.symmetry.unitcell import UnitCell


class MagneticStructure(SPWSerialisable):
    """Base class for representations of the Magnetic Structures"""

    serialisation_name = "structure"

    def __init__(self, symmetry: SymmetrySettings, super_cell: Supercell, sites: list[LatticeSite]):
        pass

    def _serialise(self, context: SPWSerialisationContext):
        return {}

    @staticmethod
    def _deserialise(data: dict, context: SPWDeserialisationContext):
        return MagneticStructure()

