from abc import ABC, abstractproperty, abstractmethod
from dataclasses import dataclass

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

class CrystalSystem(ABC):

    name: str = ""
    letter: str = ""

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

class Triclinic(CrystalSystem):
    name = "Triclinic"
    letter = "a"

    @property
    def free_parameters(self) -> FreeParameters:
        return FreeParameters()

    def constrain(self, unit_cell: UnitCell):
        return unit_cell

    @property
    def bravais_options(self) -> BravaisOptions:
        return BravaisOptions()


class Monoclinic(CrystalSystem):
    name = "Monoclinic"
    letter = "m"

    @property
    def free_parameters(self) -> FreeParameters:
        return FreeParameters(alpha=False, gamma=False)

    def constrain(self, unit_cell: UnitCell):
        return unit_cell.updated(alpha=90, beta=90)

    @property
    def bravais_options(self) -> BravaisOptions:
        return BravaisOptions(base_centered=True)


class Orthorhombic(CrystalSystem):
    name = "Orthorhombic"
    letter = "o"

    @property
    def free_parameters(self) -> FreeParameters:
        return FreeParameters(alpha=False, beta=False, gamma=False)

    def constrain(self, unit_cell: UnitCell):
        return unit_cell.updated(alpha=90, beta=90, gamma=90)

    @property
    def bravais_options(self) -> BravaisOptions:
        return BravaisOptions(
            base_centered=True,
            body_centered=True,
            face_centered=True)


class Tetragonal(CrystalSystem):
    name = "Tetragonal"
    letter = "t"

    @property
    def free_parameters(self) -> FreeParameters:
        return FreeParameters(b=False, alpha=False, beta=False, gamma=False)

    def constrain(self, unit_cell: UnitCell):
        return unit_cell.updated(
            b = unit_cell.a,
            alpha=90,
            beta=90,
            gamma=90)

    @property
    def bravais_options(self) -> BravaisOptions:
        return BravaisOptions(body_centered=True)


class Rhombohedral(CrystalSystem):
    name = "Rhombohedral"
    letter = "h"

    @property
    def free_parameters(self) -> FreeParameters:
        return FreeParameters(b=False, c=False, beta=False, gamma=False)

    def constrain(self, unit_cell: UnitCell):
        return unit_cell.updated(
            b = unit_cell.a,
            c = unit_cell.a,
            beta=unit_cell.alpha,
            gamma=unit_cell.alpha)

    @property
    def bravais_options(self) -> BravaisOptions:
        return BravaisOptions(primitive=False, rhombohedral=True)


class Hexagonal(CrystalSystem):
    name = "Hexagonal"
    letter = "h"

    @property
    def free_parameters(self) -> FreeParameters:
        return FreeParameters(b=False, alpha=False, beta=False, gamma=False)

    def constrain(self, unit_cell: UnitCell):
        return unit_cell.updated(
            b = unit_cell.a,
            alpha=90,
            beta=90,
            gamma=120)

    @property
    def bravais_options(self) -> BravaisOptions:
        return BravaisOptions()


class Cubic(CrystalSystem):
    name = "Cubic"
    letter = "c"

    @property
    def free_parameters(self) -> FreeParameters:
        return FreeParameters(b=False, c=False, alpha=False, beta=False, gamma=False)

    def constrain(self, unit_cell: UnitCell):
        return unit_cell.updated(
            b = unit_cell.a,
            c = unit_cell.a,
            alpha=90,
            beta=90,
            gamma=90)

    @property
    def bravais_options(self) -> BravaisOptions:
        return BravaisOptions(body_centered=True, face_centered=True)

# All the possible systems
crystal_systems: list[CrystalSystem] = [
    Triclinic(),
    Monoclinic(),
    Orthorhombic(),
    Tetragonal(),
    Rhombohedral(),
    Hexagonal(),
    Cubic()]

# For looking things up by name
crystal_system_name_lookup = {crystal_system.name: crystal_system for crystal_system in crystal_systems}

if __name__ == "__main__":
    for crystal_system in crystal_systems:
        symbols = [crystal_system.letter + bravais.letter for bravais in crystal_system.bravais_options.bravias]
        print(f"{crystal_system.name}: " + ",".join(symbols))

