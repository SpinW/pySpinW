""" The two types of description of the magnetic structure"""
import math
from fractions import Fraction
from typing import Iterable, Sequence

import numpy as np
from ase import Atom
from ase.lattice import BravaisLattice

from pyspinw._base import MagneticStructure
from pyspinw.checks import check_sizes


# pylint: disable=R0903


class LatticeIncomensurable(Exception):
    """ Raised when the geometry of an operation doesn't work under the constraints of a given lattice """


class MagneticLattice(MagneticStructure):
    """ A magnetic structure commensurate with the underlying crystal structure """

    def __init__(self,
                 lattice: BravaisLattice,
                 atoms: Sequence[Atom],
                 supercell_size: tuple[int, int, int] = (1,1,1)):

        super().__init__()

        self._lattice = lattice
        self._atoms = atoms
        self._supercell_size = supercell_size


    @check_sizes(vector=(3,))
    def add_propagation(self, vector: np.ndarray) -> "RotatingMagneticStructure":
        """ Create a magnetic structure based on this with an additional propagation vector"""

        return RotatingMagneticStructure(
            self._lattice,
            self._atoms,
            self._supercell_size,
            vector)

class RotatingMagneticStructure(MagneticStructure):
    """ A magnetic structure specified by a basic lattice and propagation vector

     :param lattice: Bravais lattice describing the unit cell
     :param atoms: list of atoms within the magnetic cell, with coordinates in terms of the unit cell. You might,
                   for example, have a magnetic cell made of two crystal unit cells, you would then specify
                   atoms for each crystal cell, with say, [1/2,0,0], and [3/2,0,0] on top of which you would define the
                   corresponding magnetism
     :param supercell_size: number of crystal unit cells making up the magnetic cell
     """

    @check_sizes(propagation_vector=(3,))
    def __init__(self,
                 lattice: BravaisLattice,
                 atoms: Sequence[Atom],
                 supercell_size: tuple[int, int, int] = (1,1,1),
                 propagation_vector: np.ndarray=np.array([0.0, 0.0, 0.0])):

        super().__init__()

        self._lattice = lattice
        self._atoms = atoms
        self._supercell_size = supercell_size
        self._propagation_vector = propagation_vector

    @staticmethod
    def _calculate_supercell_scaling(propagation_component: complex, supercell_component: int, max_denominator: int):
        real_frac = Fraction(propagation_component.real).limit_denominator(max_denominator)
        imag_frac = Fraction(propagation_component.imag).limit_denominator(max_denominator)

        # Cross multiply

        real = real_frac.numerator * imag_frac.denominator
        imag = imag_frac.numerator * real_frac.denominator
        cell = supercell_component * real_frac.denominator * imag_frac.denominator

        return math.lcm(real, imag, cell)

    def to_lattice(self, max_denominator: int = 10_000) -> MagneticLattice:
        """ Convert to a lattice, if possible """

        # scale each component
        new_cell_size = (RotatingMagneticStructure._calculate_supercell_scaling(
                            complex(self._propagation_vector[0]),
                            self._supercell_size[0],
                            max_denominator=max_denominator),
                         RotatingMagneticStructure._calculate_supercell_scaling(
                             complex(self._propagation_vector[1]),
                             self._supercell_size[1],
                             max_denominator=max_denominator),
                         RotatingMagneticStructure._calculate_supercell_scaling(
                             complex(self._propagation_vector[2]),
                             self._supercell_size[2],
                             max_denominator=max_denominator))