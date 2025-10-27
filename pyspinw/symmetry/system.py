""" Lattice systems

Classes and data for working with lattice systems
"""

from abc import ABC, abstractproperty, abstractmethod
from dataclasses import dataclass
from typing import Callable

from pyspinw.symmetry.bravais import PRIMITIVE, BASE_CENTERED, BODY_CENTERED, FACE_CENTERED, RHOMBOHEDRAL
from pyspinw.symmetry.unitcell import UnitCell

class IncompatibleUnitCell(Exception):
    """ Exception raised when unit cell does not fit the constraints of the lattice system"""


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

    # Note on constraints: the positive constraints describe what is needed
    #  for the unit cell to count as the given type, the negative constraints
    #  describe what is needed to not count as a different type instead.

    positive_constraints: dict[str, Callable[[UnitCell], float] | float] = {}
    """ Values of parameters that should be fixed when using this lattice type """

    negative_constraints: dict[str, Callable[[UnitCell], bool]] = {}
    """ Functions that check whether the unit cell parameters are appropriate """


    @property
    @abstractmethod
    def free_parameters(self) -> FreeParameters:
        """ List the parameters that should be available to edit"""

    def constrain(self, unit_cell: UnitCell) -> UnitCell:
        """ Constrain the unit cell for this system """
        evaluated = {
            key: value if isinstance(value, (int,float)) else value(unit_cell)
                for key, value in self.positive_constraints.items()
            }

        return unit_cell.updated(*evaluated)

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

    def validate(self, unit_cell: UnitCell):
        """ Raises an exception if the unit cell is not compatible with this system """

        # Positive constraints

        for key, value in self.positive_constraints.items():

            if not isinstance(value, (int,float)):
                value = value(unit_cell)

            if unit_cell.__dict__[key] != value:
                raise IncompatibleUnitCell(
                    f"Unit cell is not consistent with the lattice type, "
                    f"expected {key}={value}")


        # Negative constraint violations

        violations = self.violated_negative_constraints(unit_cell)
        if violations:

            # Some grammar
            if len(violations) >= 2:
                violation_string = ", ".join(violations[:-1]) + " and " + violations[-1]
            else:
                violation_string = violations[0]

            raise IncompatibleUnitCell(f"Unit cell is not consistent with the lattice type, "
                                       f"it is probably more regular than the lattice type requires. "
                                       f"Expected {violation_string}")




class Triclinic(LatticeSystem):
    """ Triclinic lattice system """

    name = "Triclinic"
    letter = "a"

    # No positive constraints

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

    @property
    def bravais_options(self) -> BravaisOptions:
        """ Bravais lattices consistent with this kind of cell"""
        return BravaisOptions()



class Monoclinic(LatticeSystem):
    """ Monoclinic lattice system """

    name = "Monoclinic"
    letter = "m"

    positive_constraints = {
        "alpha": 90,
        "beta": 90
    }

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

    @property
    def bravais_options(self) -> BravaisOptions:
        """ Bravais lattices consistent with this kind of cell"""
        return BravaisOptions(base_centered=True)

    def create_unit_cell(self, a: float, b: float, c: float, gamma: float):
        """ Create a unit cell for this lattice """

        return UnitCell(a, b, c, 90, 90, gamma)



class Orthorhombic(LatticeSystem):
    """ Orthorhombic lattice system """

    name = "Orthorhombic"
    letter = "o"

    positive_constraints = {
        "alpha": 90,
        "beta": 90,
        "gamma": 90
    }

    negative_constraints = {
        "a ≠ b": lambda cell: cell.a != cell.b,
        "b ≠ c": lambda cell: cell.b != cell.c,
        "c ≠ a": lambda cell: cell.c != cell.a,
    }

    @property
    def free_parameters(self) -> FreeParameters:
        """ List the parameters that should be available to edit"""
        return FreeParameters(alpha=False, beta=False, gamma=False)

    @property
    def bravais_options(self) -> BravaisOptions:
        """ Bravais lattices consistent with this kind of cell"""
        return BravaisOptions(
            base_centered=True,
            body_centered=True,
            face_centered=True)

    def create_unit_cell(self, a: float, b: float, c: float, gamma: float):
        """ Create a unit cell for this lattice """

        return UnitCell(a, b, c, 90, 90, gamma)

class Tetragonal(LatticeSystem):
    """ Tetragonal lattice system"""

    name = "Tetragonal"
    letter = "t"

    positive_constraints = {
        "b": lambda cell: cell.a,
        "alpha": 90,
        "beta": 90,
        "gamma": 90
    }

    negative_constraints = {
        "b ≠ c": lambda cell: cell.b != cell.c,
    }

    @property
    def free_parameters(self) -> FreeParameters:
        """ List the parameters that should be available to edit"""
        return FreeParameters(b=False, alpha=False, beta=False, gamma=False)

    @property
    def bravais_options(self) -> BravaisOptions:
        """ Bravais lattices consistent with this kind of cell"""
        return BravaisOptions(body_centered=True)


class Rhombohedral(LatticeSystem):
    """ Rhombohedral lattice system """

    name = "Rhombohedral"
    letter = "h"

    positive_constraints = {
        "b": lambda cell: cell.a,
        "c": lambda cell: cell.a,
        "beta": lambda cell: cell.alpha,
        "gamma": lambda cell: cell.alpha
    }

    negative_constraints = {
        "α ≠ 90°": lambda cell: cell.alpha != 90,
        "α ≤ 120°": lambda cell: cell.alpha <= 120
    }

    @property
    def free_parameters(self) -> FreeParameters:
        """ List the parameters that should be available to edit"""
        return FreeParameters(b=False, c=False, beta=False, gamma=False)

    @property
    def bravais_options(self) -> BravaisOptions:
        """ Bravais lattices consistent with this kind of cell"""
        return BravaisOptions(primitive=False, rhombohedral=True)


class Hexagonal(LatticeSystem):
    """ Hexagonal lattice system """

    name = "Hexagonal"
    letter = "h"

    positive_constraints = {
        "b": lambda cell: cell.a,
        "alpha": 90,
        "beta": 90,
        "gamma": 120
    }

    negative_constraints = {
        "a ≠ c": lambda cell: cell.a != cell.c,
    }

    @property
    def free_parameters(self) -> FreeParameters:
        """ List the parameters that should be available to edit"""
        return FreeParameters(b=False, alpha=False, beta=False, gamma=False)

    @property
    def bravais_options(self) -> BravaisOptions:
        """ Bravais lattices consistent with this kind of cell"""
        return BravaisOptions()


class Cubic(LatticeSystem):
    """ Cubic lattice system"""

    name = "Cubic"
    letter = "c"

    positive_constraints = {
        "b": lambda cell: cell.a,
        "c": lambda cell: cell.a,
        "alpha": 90,
        "beta": 90,
        "gamma": 90
    }

    # No negative constraints

    @property
    def free_parameters(self) -> FreeParameters:
        """ List the parameters that should be available to edit"""
        return FreeParameters(b=False, c=False, alpha=False, beta=False, gamma=False)

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

