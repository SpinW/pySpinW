""" Magnetic structure """
from pyspinw.serialisation import SPWSerialisable, SPWSerialisationContext, SPWDeserialisationContext, expects_keys
from pyspinw.site import LatticeSite
from pyspinw.symmetry.supercell import Supercell
from pyspinw.symmetry.group import SpaceGroup, MagneticSpaceGroup
from pyspinw.symmetry.unitcell import UnitCell


class MagneticStructure(SPWSerialisable):
    """Base class for representations of the Magnetic Structures"""

    serialisation_name = "structure"
    """ Object to hold symmetry information together """

    def __init__(self,
                 spacegroup: SpaceGroup,
                 magnetic_group: MagneticSpaceGroup,
                 unit_cell: UnitCell,
                 supercell: Supercell):

        self._spacegroup = spacegroup
        self._magnetic_group = magnetic_group
        self._unit_cell = unit_cell
        self._supercell = supercell

    @property
    def spacegroup(self):
        """ Get the space group """
        return self._spacegroup

    @property
    def magnetic_group(self):
        """ Get the magnetic group """
        return self._magnetic_group

    @property
    def unit_cell(self):
        """ Get the unit cell """
        return self._unit_cell

    def _serialise(self, context: SPWSerialisationContext):
        """ Serialise this object """
        return {
            "spacegroup": self._spacegroup._serialise(context),
            "magnetic_group": self._magnetic_group._serialise(context),
            "unit_cell": self._unit_cell._serialise(context),
            "supercell": self._supercell._serialise(context)
        }

    @staticmethod
    @expects_keys("spacegroup, magnetic_group, unit_cell, supercell")
    def _deserialise(data: dict, context: SPWDeserialisationContext):
        """ Deserialise an object of this kind """
        return MagneticStructure(
            spacegroup=SpaceGroup._deserialse(data["spacegroup"]),
            magnetic_group=MagneticSpaceGroup._deserialise(data["magnetic_group"]),
            unit_cell=UnitCell._deserialise(data["unit_cell"]),
            supercell=Supercell._deserialise(data["supercell"])
        )

    def validate(self):
        """ Check that sites are compatible with the symmetry groups and supercell """
        # TODO: Implement
        raise NotImplementedError()
