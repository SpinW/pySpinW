""" Classes for working with spacegroup settings """

import re
from dataclasses import dataclass
from enum import Enum


class RhombohedralOrHexagonal(Enum):
    """ Rhombohedral or hexagonal choice """

    RHOMBOHEDRAL = "R"
    HEXAGONAL = "H"

class CrystalAxis(Enum):
    """ A crystal axis"""

    A = "a"
    B = "b"
    C = "c"

@dataclass
class UniqueAxis:
    """ Description of the unique axis of spacegroup"""

    axis: CrystalAxis
    negated: bool = False

    def to_string(self):
        """ Convert to spglib type string """
        s = "-" if self.negated else ""
        s += self.axis.value
        return s

    @staticmethod
    def from_string(string: str):
        """ Create from spglib type string"""
        negated = string.startswith("-")

        axis_string = string[1] if negated else string[0]

        axis = CrystalAxis(axis_string)

        return UniqueAxis(axis, negated)



class AxisPermutation:
    """ Description of a permutation of a unit cell """

    def __init__(self,
                 first_axis: CrystalAxis,
                 second_axis: CrystalAxis,
                 flip_c: bool = False):

        if first_axis == second_axis:
            raise ValueError("Expected first and second axes to be different")

        self.first_axis = first_axis
        self.second_axis = second_axis

        third = [CrystalAxis.A, CrystalAxis.B, CrystalAxis.C]
        third.remove(first_axis)
        third.remove(second_axis)

        self.third_axis = third[0]

        self.flip_c = flip_c

    def to_string(self):
        """ Convert permutation to string """
        s = ""
        for axis in [self.first_axis, self.second_axis, self.third_axis]:
            if axis == CrystalAxis.C and self.flip_c:
                s += "-"
            s += axis.value

        return s

    @staticmethod
    def from_string(string: str):
        """ Convert a string e.g. "b-ca" to a permutation representation """
        axes = []
        flip = False
        # Parse the string
        for char in string:
            match char:
                case "a":
                    axes.append(CrystalAxis.A)
                case "b":
                    axes.append(CrystalAxis.B)
                case "c":
                    axes.append(CrystalAxis.C)
                case "-":
                    flip = True
                case _:
                    raise ValueError("Expected string to contain only the characters 'a', 'b', 'c' and '-'")

        perm = AxisPermutation(axes[0], axes[1], flip_c=flip)

        if perm.to_string() != string:
            raise ValueError("Malformed input string")

        return perm

@dataclass
class Setting:
    """ Spacegroup setting descriptor in a form that can be easily worked with """

    _re_pattern = (
        r'^(?P<rhombhex>[RH])?'  # Rhombohedral or Hexagonal
        r'(?P<unique>-?[abc])?'  # Unique Axis
        r'(?P<choice>\d)?'  # Choice (single digit)
        r'(?P<perm>(?:-?c|[ab]){3})?$'  # Permutation
    )

    choice_number: int | None = None,
    permutation: AxisPermutation | None = None,
    unique_axis: UniqueAxis | None = None,
    rhombohedral_or_hexagonal: RhombohedralOrHexagonal | None = None

    @staticmethod
    def from_optional_string(string: str | None):
        """ Convert from spglib database string, or default setting if None"""
        if string is None:
            return Setting()
        else:
            return Setting.from_string(string)

    @staticmethod
    def from_string(string: str):
        """ Convert from spglib database string"""
        m = re.match(Setting._re_pattern, string)

        if not m:
            raise ValueError("Not a valid setting string")

        rhombhex = None if m["rhombhex"] is None else RhombohedralOrHexagonal(m["rhombhex"])
        unique = None if m["unique"] is None else UniqueAxis.from_string(m["unique"])
        choice = None if m["choice"] is None else int(m["choice"])
        permutation = None if m["perm"] is None else AxisPermutation.from_string(m["perm"])

        return Setting(choice, permutation, unique, rhombohedral_or_hexagonal=rhombhex)

    def to_string(self):
        """ Convert this setting to a string like it is in the spglib database """
        s = "" if self.rhombohedral_or_hexagonal is None else self.rhombohedral_or_hexagonal.value
        s += "" if self.unique_axis is None else self.unique_axis.to_string()
        s += "" if self.choice_number is None else str(self.choice_number)
        s += "" if self.permutation is None else self.permutation.to_string()
        return s
