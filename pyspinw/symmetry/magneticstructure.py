
from pyspinw.serialisation import SPWSerialisable, SPWSerialisationContext, SPWDeserialisationContext
from pyspinw.symmetry.supercell import Supercell
from pyspinw.symmetry.group import SpaceGroup, MagneticSpaceGroup
from pyspinw.symmetry.unitcell import UnitCell


class MagneticStructure(SPWSerialisable):
    """Base class for representations of the Magnetic Structures"""

    serialisation_name = "structure"
    """ Object to hold symmetry information together """

    def __init__(self, space_group: SpaceGroup, magnetic_group: MagneticSpaceGroup, unit_cell: UnitCell, supercell: Supercell):
        self._space_group = space_group
        self._magnetic_group = magnetic_group
        self._unit_cell = unit_cell
        self._supercell = supercell

    @property
    def space_group(self):
        return self._space_group

    @property
    def magnetic_group(self):
        return self._magnetic_group

    @property
    def unit_cell(self):
        return self._unit_cell

    def _serialise(self, context: SPWSerialisationContext):
        return {}

    @staticmethod
    def _deserialise(data: dict, context: SPWDeserialisationContext):
        return MagneticStructure()

    def validate(self):
        """ Check that sites are compatible with the symmetry groups and supercell """
        # TODO: Implement
        raise NotImplementedError()
