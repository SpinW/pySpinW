""" Lattice systems

Classes and data for working with lattice systems
"""

from abc import ABC, abstractproperty, abstractmethod
from dataclasses import dataclass
from typing import Callable

from pyspinw.symmetry.bravais import PRIMITIVE, BASE_CENTERED, BODY_CENTERED, FACE_CENTERED, RHOMBOHEDRAL
from pyspinw.symmetry.unitcell import UnitCell


@dataclass
class FreeParameters:
    """ Which parameters of the unit cell are free to vary """

    a: bool = True
    b: bool = True
    c: bool = True
    alpha: bool = True
    beta: bool = True
    gamma: bool = True


@dataclass
class BravaisOptions:
    """ Types of Bravias lattices available """

    primitive: bool = True
    base_centered: bool = False
    body_centered: bool = False
    face_centered: bool = False
    rhombohedral: bool = False

    @property
    def bravias(self):
        """List all possible uppercase lattice letters"""
        out = []

        if self.primitive:
            out.append(PRIMITIVE)

        if self.base_centered:
            out.append(BASE_CENTERED)

        if self.body_centered:
            out.append(BODY_CENTERED)

        if self.face_centered:
            out.append(FACE_CENTERED)

        if self.rhombohedral:
            out.append(RHOMBOHEDRAL)

        return out

class LatticeSystem(ABC):
    """ Base class for the various lattice systems """

    name: str = ""
    letter: str = ""

    negative_constraints: dict[str, Callable[[UnitCell], bool]] = {}
    """ Functions that check whether the unit cell parameters are appropriate """

    @property
    @abstractmethod
    def free_parameters(self) -> FreeParameters:
        """ List the parameters that should be available to edit"""

    @abstractmethod
    def constrain(self, unit_cell: UnitCell) -> UnitCell:
        """ Constrain the unit cell for this system """

    @property
    @abstractmethod
    def bravais_options(self) -> BravaisOptions:
        """ Bravais lattices consistent with this kind of cell"""

    def violated_negative_constraints(self, unit_cell: UnitCell) -> list[str]:
        """ Check whether any of the constraints on specific values are violated by this unit cell """
        out = []
        for name, constraint in self.negative_constraints.items():
            if not constraint(unit_cell):
                out.append(name)
        return out


class Triclinic(LatticeSystem):
    """ Triclinic lattice system """

    name = "Triclinic"
    letter = "a"

    negative_constraints = {
        "a ≠ b": lambda cell: cell.a != cell.b,
        "b ≠ c": lambda cell: cell.b != cell.c,
        "c ≠ a": lambda cell: cell.c != cell.a,
        "α ≠ 90°": lambda cell: cell.alpha != 90,
        "β ≠ 90°": lambda cell: cell.beta != 90,
        "γ ≠ 90°": lambda cell: cell.gamma != 90
    }

    @property
    def free_parameters(self) -> FreeParameters:
        """ List the parameters that should be available to edit"""
        return FreeParameters()

    def constrain(self, unit_cell: UnitCell):
        """ Constrain the unit cell for this system """
        return unit_cell

    @property
    def bravais_options(self) -> BravaisOptions:
        """ Bravais lattices consistent with this kind of cell"""
        return BravaisOptions()



class Monoclinic(LatticeSystem):
    """ Monoclinic lattice system """

    name = "Monoclinic"
    letter = "m"

    negative_constraints = {
        "a ≠ b": lambda cell: cell.a != cell.b,
        "b ≠ c": lambda cell: cell.b != cell.c,
        "c ≠ a": lambda cell: cell.c != cell.a,
        "γ ≠ 90°": lambda cell: cell.gamma != 90
    }

    @property
    def free_parameters(self) -> FreeParameters:
        """ List the parameters that should be available to edit"""
        return FreeParameters(alpha=False, gamma=False)

    def constrain(self, unit_cell: UnitCell):
        """ Constrain the unit cell for this system """
        return unit_cell.updated(alpha=90, beta=90)

    @property
    def bravais_options(self) -> BravaisOptions:
        """ Bravais lattices consistent with this kind of cell"""
        return BravaisOptions(base_centered=True)


class Orthorhombic(LatticeSystem):
    """ Orthorhombic lattice system """

    name = "Orthorhombic"
    letter = "o"


    negative_constraints = {
        "a ≠ b": lambda cell: cell.a != cell.b,
        "b ≠ c": lambda cell: cell.b != cell.c,
        "c ≠ a": lambda cell: cell.c != cell.a,
    }

    @property
    def free_parameters(self) -> FreeParameters:
        """ List the parameters that should be available to edit"""
        return FreeParameters(alpha=False, beta=False, gamma=False)

    def constrain(self, unit_cell: UnitCell):
        """ Constrain the unit cell for this system """
        return unit_cell.updated(alpha=90, beta=90, gamma=90)

    @property
    def bravais_options(self) -> BravaisOptions:
        """ Bravais lattices consistent with this kind of cell"""
        return BravaisOptions(
            base_centered=True,
            body_centered=True,
            face_centered=True)


class Tetragonal(LatticeSystem):
    """ Tetragonal lattice system"""

    name = "Tetragonal"
    letter = "t"

    negative_constraints = {
        "b ≠ c": lambda cell: cell.b != cell.c,
    }

    @property
    def free_parameters(self) -> FreeParameters:
        """ List the parameters that should be available to edit"""
        return FreeParameters(b=False, alpha=False, beta=False, gamma=False)

    def constrain(self, unit_cell: UnitCell):
        """ Constrain the unit cell for this system """
        return unit_cell.updated(
            b = unit_cell.a,
            alpha=90,
            beta=90,
            gamma=90)

    @property
    def bravais_options(self) -> BravaisOptions:
        """ Bravais lattices consistent with this kind of cell"""
        return BravaisOptions(body_centered=True)


class Rhombohedral(LatticeSystem):
    """ Rhombohedral lattice system """

    name = "Rhombohedral"
    letter = "h"

    negative_constraints = {
        "α ≠ 90°": lambda cell: cell.alpha != 90,
        "α ≤ 120°": lambda cell: cell.alpha <= 120
    }

    @property
    def free_parameters(self) -> FreeParameters:
        """ List the parameters that should be available to edit"""
        return FreeParameters(b=False, c=False, beta=False, gamma=False)

    def constrain(self, unit_cell: UnitCell):
        """ Constrain the unit cell for this system """
        return unit_cell.updated(
            b = unit_cell.a,
            c = unit_cell.a,
            beta=unit_cell.alpha,
            gamma=unit_cell.alpha)

    @property
    def bravais_options(self) -> BravaisOptions:
        """ Bravais lattices consistent with this kind of cell"""
        return BravaisOptions(primitive=False, rhombohedral=True)


class Hexagonal(LatticeSystem):
    """ Hexagonal lattice system """

    name = "Hexagonal"
    letter = "h"

    negative_constraints = {
        "a ≠ c": lambda cell: cell.a != cell.c,
    }

    @property
    def free_parameters(self) -> FreeParameters:
        """ List the parameters that should be available to edit"""
        return FreeParameters(b=False, alpha=False, beta=False, gamma=False)

    def constrain(self, unit_cell: UnitCell):
        """ Constrain the unit cell for this system """
        return unit_cell.updated(
            b = unit_cell.a,
            alpha=90,
            beta=90,
            gamma=120)

    @property
    def bravais_options(self) -> BravaisOptions:
        """ Bravais lattices consistent with this kind of cell"""
        return BravaisOptions()


class Cubic(LatticeSystem):
    """ Cubic lattice system"""

    name = "Cubic"
    letter = "c"

    @property
    def free_parameters(self) -> FreeParameters:
        """ List the parameters that should be available to edit"""
        return FreeParameters(b=False, c=False, alpha=False, beta=False, gamma=False)

    def constrain(self, unit_cell: UnitCell):
        """ Constrain the unit cell for this system """
        return unit_cell.updated(
            b = unit_cell.a,
            c = unit_cell.a,
            alpha=90,
            beta=90,
            gamma=90)

    @property
    def bravais_options(self) -> BravaisOptions:
        """ Bravais lattices consistent with this kind of cell"""
        return BravaisOptions(body_centered=True, face_centered=True)

# All the possible systems
lattice_systems: list[LatticeSystem] = [
    Triclinic(),
    Monoclinic(),
    Orthorhombic(),
    Tetragonal(),
    Rhombohedral(),
    Hexagonal(),
    Cubic()]

# For looking things up by name or letter
lattice_system_name_lookup = {lattice.name: lattice for lattice in lattice_systems}
lattice_system_letter_lookup = {lattice.letter: lattice for lattice in lattice_systems}

def find_unit_cell_type(unit_cell: UnitCell) -> list[LatticeSystem]:
    """ Find the lattice system for a unit cell, it should be a length zero or one list, but it might not be """
    # get potential lattice systems
    potential_systems = [system for system in lattice_systems if system.constrain(unit_cell) == unit_cell]

    # check for violated constraints, and return
    return [system for system in potential_systems if len(system.violated_negative_constraints(unit_cell)) == 0]




if __name__ == "__main__":
    for lattice_system in lattice_systems:
        symbols = [lattice_system.letter + bravais.letter for bravais in lattice_system.bravais_options.bravias]
        print(f"{lattice_system.name}: " + ",".join(symbols))

