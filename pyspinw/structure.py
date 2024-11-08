""" The two types of description of the magnetic structure"""

from pyspinw._base import MagneticStructure


class LatticeIncomensurable(Exception):
    """ Raised when the geometry of an operation doesn't work under the constraints of a given lattice """


class MagneticLattice(MagneticStructure):
    """ A magnetic structure commensurate with the underlying crystal structure """


    def add_propagation(self) -> "MagneticPsuedolattice":
        """ Create a magnetic structure based on this with an additional propagation vector"""


class MagneticPsuedolattice(MagneticStructure):
    """ A magnetic structure specified by a basic lattice and propagation vector """

    def to_lattice(self) -> MagneticLattice:
        """ Convert to a lattice, if possible """