""" Things needed specifically for deserialisation, basically a list of objects that can be serialised """
from pyspinw import SiteMetadata, LatticeSite, UnitCell, Anisotropy, Structure
from pyspinw.exchange import Exchange
from pyspinw.sample import Sample
from pyspinw.site import ImpliedLatticeSite
from pyspinw.symmetry.group import SpaceGroup
from pyspinw.symmetry.supercell import Supercell

# A list of the classes that can be specified in the serialised json,
# this is needed so we can load objects of different kinds
serialisation_entry_classes = [
    SiteMetadata,
    LatticeSite,
    ImpliedLatticeSite,
    SpaceGroup,
    UnitCell,
    Structure,
    Exchange,
    Anisotropy,
    Supercell,
    Sample
]