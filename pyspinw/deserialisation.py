"""Things needed specifically for deserialisation, basically a list of objects that can be serialised.

Attributes
----------
serialisation_entry_classes : list of type
    List of classes that can be specified in serialised JSON.
serialisation_class_lookup : dict of {str: type}
    Mapping from ``serialisation_name`` string to serialisable class.
"""
from pyspinw import CrystalDomain
from pyspinw.cell_offsets import CellOffset
from pyspinw.exchange import Exchange
from pyspinw.exchangemetadata import ExchangeMetadata
from pyspinw.experiment import Experiment
from pyspinw.sample import Sample
from pyspinw.site import SiteMetadata, Site, ImpliedLatticeSite
from pyspinw.symmetry.group import SpaceGroup
from pyspinw.symmetry.magnetic_group import MagneticSpaceGroup
from pyspinw.symmetry.supercell import Supercell, PropagationVector, SupercellTransformation
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structure import Structure
from pyspinw.anisotropy import Anisotropy
from pyspinw.hamiltonian import Hamiltonian

# A list of the classes that can be specified in the serialised json,
# this is needed so we can load objects of different kinds
serialisation_entry_classes = [
    SiteMetadata,
    Site,
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

    CrystalDomain,
    Sample,

    Experiment
]

serialisation_class_lookup = {cls.serialisation_name: cls for cls in serialisation_entry_classes}
