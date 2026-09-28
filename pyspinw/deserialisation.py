""" Things needed specifically for deserialisation, basically a list of objects that can be serialised """

from pyspinw.cell_offsets import CellOffset
from pyspinw.exchange import Exchange
from pyspinw.exchangemetadata import ExchangeMetadata
from pyspinw.experiment import Experiment
from pyspinw.sample import Sample
from pyspinw.site import SiteMetadata, LatticeSite, ImpliedLatticeSite
from pyspinw.symmetry.group import SpaceGroup, MagneticSpaceGroup
from pyspinw.symmetry.supercell import Supercell, PropagationVector, SupercellTransformation
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structure import Structure
from pyspinw.anisotropy import Anisotropy
from pyspinw.hamiltonian import Hamiltonian

# A list of the classes that can be specified in the serialised json,
# this is needed so we can load objects of different kinds
serialisation_entry_classes = [
    SiteMetadata,
    LatticeSite,
    ImpliedLatticeSite,

    SpaceGroup,
    MagneticSpaceGroup,
    UnitCell,
    SupercellTransformation,
    PropagationVector,
    Supercell,
    Structure,

    CellOffset,
    ExchangeMetadata,
    Exchange,
    Anisotropy,

    Hamiltonian,

    Sample,

    Experiment
]

serialisation_class_lookup = {cls.serialisation_name: cls for cls in serialisation_entry_classes}