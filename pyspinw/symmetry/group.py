""" Space groups and magnetic space groups"""

from collections import defaultdict
from typing import Callable

import numpy as np
import spglib
from difflib import get_close_matches

from pyspinw.site import LatticeSite, ImpliedLatticeSite
from pyspinw.symmetry.group_conventions import spacegroup_conventions, canonise_spacegroup_name

from pyspinw.symmetry.operations import MagneticOperation, SpaceOperation
from pyspinw.symmetry.data.msg_symbols import msg_symbols
from pyspinw.symmetry.system import LatticeSystem, lattice_system_letter_lookup, Rhombohedral, \
    lattice_system_name_lookup

from pyspinw.tolerances import tolerances


class SymmetryGroup:
    """Base class for the different kinds of groups """

class MagneticSpaceGroup(SymmetryGroup):
    """ Representation of a magnetic space group"""

    def __init__(self, number: int, symbol: str, operations: list[MagneticOperation]):
        self.number = number
        self.symbol = symbol
        self.operations = operations

    def __repr__(self):
        """repr"""
        return f"MagneticSpaceGroup({self.number}, {self.symbol})"

    def duplicates(self, site: LatticeSite) -> list[ImpliedLatticeSite]:
        """ Find "duplicate" sites of a given site """
        coordinates = site.values.reshape(1, -1)

        new_coordinates = []
        for operation in self.operations:

            candidate = operation(coordinates)

            # If its not the input, continue
            if np.all(np.abs(candidate - coordinates) < tolerances.SAME_SITE_ABS_TOL):
                continue

            # Is it one we've already found
            new = True
            for ijkm in new_coordinates:
                if np.all(np.abs(candidate - ijkm) < tolerances.SAME_SITE_ABS_TOL):
                    new = False
                    break

            if new:
                new_coordinates.append(candidate)


        new_sites = []
        for i, ijkm in enumerate(new_coordinates):
            new_site = ImpliedLatticeSite.from_coordinates(
                coordinates=ijkm.reshape(-1),
                parent_site=site,
                name=site.name + f" [{i+1}]")

            new_sites.append(new_site)

        return new_sites



class SpaceGroup[T: LatticeSystem](SymmetryGroup):
    """ Representation of a space group"""

    def __init__(self,
                 number,
                 international_symbol,
                 short_symbol,
                 operations,
                 magnetic_variants: list[MagneticSpaceGroup],
                 lattice_system: T,
                 choice: str | None):

        self.number = number
        self.symbol = international_symbol
        self.short_symbol = short_symbol
        self.operations = operations
        self.magnetic_variants = magnetic_variants
        self.lattice_system: T = lattice_system
        self.choice = choice

        # This is slightly unusual, make a reference to the lattice system create_unit_cell method
        self.create_unit_cell = lattice_system.create_unit_cell

    def __repr__(self):
        """repr"""
        if self.choice is None:
            return f"SpaceGroup({self.number}, {self.symbol})"
        else:
            return f"SpaceGroup({self.number}, {self.symbol} [{self.choice}])"


class NoSuchGroup(Exception):
    """ Raised when a group is not found. """



class SpacegroupDatabase:
    """ Holds all the information about spacegroups and magnetic spacegroups """

    # Letter representing the crystal system monoclinic, cubic etc
    _spacegroup_number_to_lattice_system = \
        ["a" for _ in range(1, 3)] + \
        ["m" for _ in range(3, 16)] + \
        ["o" for _ in range(16, 75)] + \
        ["t" for _ in range(75, 143)] + \
        ["h" for _ in range(143, 195)] + \
        ["c" for _ in range(195, 231)]

    # Converts first letter of spacegroup to the relevant Bravais notation
    _spacegroup_symbol_to_bravais_symbol = {
        "A": "S",
        "B": "S",
        "C": "S",
        "F": "F",
        "I": "I",
        "P": "P",
        "R": "R",
    }

    # Used for looking up spacegroups by name
    _canonincal_name_to_group_name, _canonical_spacegroup_name_to_index = spacegroup_conventions()

    def __init__(self):

        # We need to know which magnetic group is associated with each spacegroup before
        #  constructing the spacegroup instances
        spacegroup_to_magnetic_group = defaultdict(list[int])
        for i in range(1,1652):
            metadata = spglib.get_magnetic_spacegroup_type(i)
            spacegroup_to_magnetic_group[metadata.number].append(i)

        # Create magnetic groups
        self.magnetic_groups = []
        for i in range(1,1652):

            op_data = spglib.get_magnetic_symmetry_from_database(i)
            # metadata = spglib.get_magnetic_spacegroup_type(i)
            # print(metadata)

            rotations = op_data["rotations"]
            translations = op_data["translations"]
            time_reversals = 1 - 2*op_data["time_reversals"]

            symbol = msg_symbols[i].uni

            operations = []
            for rotation, translation, time_reversal in zip(rotations, translations, time_reversals):
                op = MagneticOperation.from_numpy(rotation, translation, time_reversal)
                operations.append(op)

            group = MagneticSpaceGroup(i, symbol, operations)
            self.magnetic_groups.append(group)



        # Create spacegroups
        self._lattice_symbol_to_spacegroups = defaultdict(list[SpaceGroup])
        self.spacegroups: list[SpaceGroup] = []


        for i in range(1, 531):

            # Get the relevant data

            sg_data = spglib.get_spacegroup_type(i)
            op_data = spglib.get_symmetry_from_database(i)
            corresponding_magnetic_groups = [self.magnetic_groups[idx-1] for idx in spacegroup_to_magnetic_group[i]]

            # Classify
            name = sg_data.international_full
            short_name = sg_data.international_short
            lattice_type_from_name = name[0]

            first_letter = SpacegroupDatabase._spacegroup_number_to_lattice_system[sg_data.number-1]
            second_letter = SpacegroupDatabase._spacegroup_symbol_to_bravais_symbol[lattice_type_from_name]
            bravais_lattice = first_letter + second_letter

            # Load the operations
            translations = op_data["translations"]
            rotations = op_data["rotations"]

            operations = []
            for translation, rotation, in zip(translations, rotations):
                op = SpaceOperation.from_numpy(rotation, translation)
                operations.append(op)

            choice = None if sg_data.choice == "" else sg_data.choice

            if choice == "R":
                lattice_system = lattice_system_name_lookup["Rhombohedral"]
                bravais_lattice = "hR"
            else:
                lattice_system = lattice_system_letter_lookup[first_letter]

            group = SpaceGroup(
                number=sg_data.number,
                international_symbol=name,
                short_symbol=short_name,
                operations=operations,
                choice=choice,
                magnetic_variants=corresponding_magnetic_groups,
                lattice_system=lattice_system)

            self.spacegroups.append(group)
            self._lattice_symbol_to_spacegroups[bravais_lattice].append(group)

        # Lookup for spacegroups by long/short international symbol
        self._spacegroup_symbol_lookup: dict[str, SpaceGroup] = {group.symbol: group for group in self.spacegroups}
        self._spacegroup_symbol_lookup.update({group.short_symbol: group for group in self.spacegroups})

        # Lookup for magnetic spacegroups
        self._magnetic_group_symbol_lookup: dict[str, MagneticSpaceGroup] = {
                group.symbol: group for group in self.magnetic_groups }

        self._canonical_magnetic_spacegroup_name_to_index = {
            canonise_spacegroup_name(group.symbol): group for group in self.magnetic_groups }

        # Lookup for spacegroups based on their number
        self._spacegroup_number_lookup = defaultdict(list[SpaceGroup])
        for group in self.spacegroups:
            self._spacegroup_number_lookup[group.number].append(group)



    def all_spacegroup_settings(self, number: int):
        """ Get a list of all spacegroups that correspond to a given setting """
        return self._spacegroup_number_lookup[number]

    def spacegroups_for_lattice_symbol(self, lattice_symbol: str):
        return self._lattice_symbol_to_spacegroups[lattice_symbol]

    def spacegroup_by_name(self, name: str) -> SpaceGroup:
        """ Get a spacegroup by name"""

        canonised_input = canonise_spacegroup_name(name)

        if canonised_input in self._canonical_spacegroup_name_to_index:
            index = self._canonical_spacegroup_name_to_index[canonised_input]
            return self.spacegroups[index-1]

        # Failure branch...
        # We want to get names that are similar in their input text form, but report the standard form
        similar = get_close_matches(name, self._canonincal_name_to_group_name.keys(), n=3, cutoff=0.4)

        if similar:
            suggestion_string = ", ".join([f"'{self._canonincal_name_to_group_name[s]}'" for s in similar])
            message_string = f"Unknown space group '{name}', perhaps you meant {suggestion_string} or something similar."

        else:
            message_string = f"Unknown space group '{name}', doesn't seem to be even close to a magnetic space group."

        raise NoSuchGroup(message_string)


    def magnetic_spacegroup_by_name(self, name: str) -> MagneticSpaceGroup:
        """ Get a magnetic spacegroup by name

        TODO: make this work more like the normal spacegroup one
        """
        #
        # lower_name = name.lower()
        # if lower_name in _lowercase_magnetic_group_lookup:
        #     return _lowercase_magnetic_group_lookup[lower_name]
        #
        # similar = get_close_matches(lower_name, _lowercase_magnetic_group_lookup.keys(), n=3)
        # formatted = [_lowercase_magnetic_group_lookup[key].symbol for key in similar]
        #
        # if similar:
        #     suggestion_string = ", ".join([f"'{s}'" for s in formatted])
        #     message_string = f"Unknown space group '{name}', perhaps you meant {suggestion_string} or something similar."
        #
        # else:
        #
        #     message_string = f"Unknown space group '{name}', doesn't seem to be even close to a magnetic space group."
        #
        # raise NoSuchGroup(message_string)


database = SpacegroupDatabase()

if __name__ == "__main__":
    print(database.spacegroup_by_name("b2/m"))
    print(database.spacegroup_by_name("c2/m"))

    print(database.spacegroup_by_name("pnnn:1"))

    for i, spacegroup in enumerate(database.spacegroups):
        print(i, spacegroup)
