""" Space groups and magnetic space groups"""

from collections import defaultdict
from enum import Enum
from typing import Callable

from abc import ABC, abstractmethod

import numpy as np
import spglib
from difflib import get_close_matches

from pyspinw.serialisation import SPWSerialisable, SPWSerialisationContext, SPWDeserialisationContext
from pyspinw.site import LatticeSite, ImpliedLatticeSite
from pyspinw.symmetry.canonise import canonise_string
from pyspinw.symmetry.spacegroup_lookup import canonical_aliases, canonical_to_formatted, preferred_names

from pyspinw.symmetry.operations import MagneticOperation, SpaceOperation
from pyspinw.symmetry.data.msg_symbols import msg_symbols
from pyspinw.symmetry.settings import Setting
from pyspinw.symmetry.supercell import Supercell
from pyspinw.symmetry.system import LatticeSystem, lattice_system_letter_lookup, Rhombohedral, \
    lattice_system_name_lookup

from pyspinw.tolerances import tolerances


class SymmetryGroup(ABC, SPWSerialisable):
    """ Base class for symmetry group and magnetic symmetry group """

    @abstractmethod
    def implied_sites_for(self, site: LatticeSite) -> list[ImpliedLatticeSite]:
        """ Find all the sites that are required by symmetry by the input site """




class MagneticSpaceGroup(SymmetryGroup):
    """ Representation of a magnetic space group"""

    serialisation_name = "MagneticGroup"

    def __init__(self, number: int, symbol: str, operations: list[MagneticOperation]):
        self.number = number
        self.symbol = symbol
        self.operations = operations

    def __repr__(self):
        """repr"""
        return f"MagneticSpaceGroup({self.number}, {self.symbol})"

    def implied_sites_for(self, site: LatticeSite) -> list[ImpliedLatticeSite]:
        """ Find "duplicate" sites of a given site """
        # TODO: Transform moments in the correct coordinate system

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



class SpaceGroup(SymmetryGroup):
    """ Representation of a space group"""

    serialisation_name = "SpaceGroup"

    def __init__(self,
                 number,
                 international_symbol,
                 short_symbol,
                 preferred_symbol,
                 operations,
                 magnetic_variants: list[MagneticSpaceGroup],
                 lattice_system: LatticeSystem,
                 choice: str | None):

        self.number = number
        self.symbol = international_symbol
        self.short_symbol = short_symbol
        self.preferred_symbol = preferred_symbol
        self.operations = operations
        self.magnetic_variants = magnetic_variants
        self.lattice_system = lattice_system
        self.choice = choice
        self.setting = Setting.from_optional_string(choice)

        # This is slightly unusual, make a reference to the lattice system create_unit_cell method
        self.create_unit_cell = lattice_system.create_unit_cell

    def __repr__(self):
        """repr"""
        if self.choice is None:
            return f"SpaceGroup({self.number}, {self.symbol})"
        else:
            return f"SpaceGroup({self.number}, {self.symbol} [{self.choice}])"

    def _serialisation_string(self):
        """ Name to use to refer to this group in serialisation"""
        return self.preferred_symbol

    def _serialise(self, context: SPWSerialisationContext):
        pass

    @staticmethod
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        pass

    def for_supercell(self, supercell: Supercell):
        """ Get the symmetry group of a supercell, as implied by the symmetry of the unit cell """

        return database.spacegroup_by_name("p1")

    def implied_sites_for(self, site: LatticeSite) -> list[ImpliedLatticeSite]:
        """ Find "duplicate" sites of a given site """
        coordinates = site.ijk.reshape(1, -1)

        new_coordinates = []
        for operation in self.operations:

            candidate = operation(coordinates)

            # If it's not the input, continue
            if np.all(np.abs(candidate - coordinates) < tolerances.SAME_SITE_ABS_TOL):
                continue

            # Is it one we've already found
            new = True
            for coordinates in new_coordinates:
                if np.all(np.abs(candidate - coordinates) < tolerances.SAME_SITE_ABS_TOL):
                    new = False
                    break

            if new:
                new_coordinates.append(candidate)


        new_sites = []
        for i, coordinates in enumerate(new_coordinates):

            new_site = ImpliedLatticeSite(
                parent_site=site,
                i=coordinates[0][0],
                j=coordinates[0][1],
                k=coordinates[0][2],
                supercell_moments=site._moment_data,
                name=site.name + f" [{i+1}]"
                )

            new_sites.append(new_site)

        return new_sites

#
#
# Database stuff
#
#

class ExactMatch:
    """ Result for exact matches in database searches"""

    def __init__(self, spacegroup: SpaceGroup):
        self.spacegroup = spacegroup

class PartialMatch:
    """ Result for database searches that are not exact"""

    def __init__(self, spacegroups: list[SpaceGroup]):
        self.spacegroups = spacegroups


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
        self._operation_lookup = []

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

            self._operation_lookup.append(
                {canonise_string(operation.text_form) for operation in operations})

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
                preferred_symbol=preferred_names[i],
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
            canonise_string(group.symbol): group for group in self.magnetic_groups }

        # Lookup for spacegroups based on their number
        self._spacegroup_number_lookup = defaultdict(list[SpaceGroup])
        for group in self.spacegroups:
            self._spacegroup_number_lookup[group.number].append(group)



    def all_spacegroup_settings(self, number: int):
        """ Get a list of all spacegroups that correspond to a given setting """
        return self._spacegroup_number_lookup[number]

    def spacegroups_for_lattice_symbol(self, lattice_symbol: str):
        """ Get the spacegroups that correspond to a given lattice symbol (e.g. hR)"""
        return self._lattice_symbol_to_spacegroups[lattice_symbol]

    def spacegroups_with_operations(self, operation_string: str) -> ExactMatch | PartialMatch:
        """ Get a spacegroup by a list of operations

        Example: TODO

        :param operation_string: A list of operations in the x,y,z form, semicolon separated
        :return: ExactMatch or MultipleMatches
        """
        operations = {canonise_string(op) for op in operation_string.split(";")}

        out = []
        for i, test_operations in enumerate(self._operation_lookup):
            if operations <= test_operations:
                if operations == test_operations:
                    return ExactMatch(self.spacegroups[i])

                out.append(self.spacegroups[i])

        return PartialMatch(out)


    def spacegroup_by_name(self, name: str) -> SpaceGroup:
        """ Get a spacegroup by name"""
        canonised_input = canonise_string(name)

        if canonised_input in canonical_aliases:
            index = canonical_aliases[canonised_input]
            return self.spacegroups[index-1]

        # Failure branch...
        # We want to get names that are similar in their input text form, but report the standard form
        similar = get_close_matches(name, canonical_aliases.keys(), n=3, cutoff=0.4)

        if similar:
            suggestion_string = ", ".join([f"'{canonical_to_formatted[s]}'" for s in similar])
            message_string = (f"Unknown space group '{name}', "
                              f"perhaps you meant {suggestion_string} or something similar.")

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
        #     message_string = (f"Unknown space group '{name}', "
        #                       f"perhaps you meant {suggestion_string} or something similar.")
        #
        # else:
        #
        #     message_string = f"Unknown space group '{name}', doesn't seem to be even close to a magnetic space group."
        #
        # raise NoSuchGroup(message_string)


database = SpacegroupDatabase()

if __name__ == "__main__":
    # print(database.spacegroup_by_name("b2/m"))
    print(database.spacegroup_by_name("c2/m"))

    print(database.spacegroup_by_name("pnnn:1"))

    for i, spacegroup in enumerate(database.spacegroups):
        print(i, spacegroup)
