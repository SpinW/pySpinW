""" Space groups and magnetic space groups"""

from collections import defaultdict

import numpy as np
import spglib
from difflib import get_close_matches

from pyspinw.site import LatticeSite, ImpliedLatticeSite
from pyspinw.symmetry.operations import MagneticOperation, SpaceOperation
from pyspinw.symmetry.data.msg_symbols import msg_symbols
from pyspinw.tolerances import tolerances
from symmetry.system import LatticeSystem, lattice_system_letter_lookup


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



class SpaceGroup(SymmetryGroup):
    """ Representation of a spacegroup"""

    def __init__(self,
                 number,
                 international_symbol,
                 short_symbol,
                 operations,
                 magnetic_variants: list[MagneticSpaceGroup],
                 lattice_system: LatticeSystem):

        self.number = number
        self.symbol = international_symbol
        self.short_symbol = short_symbol
        self.operations = operations
        self.magnetic_variants = magnetic_variants
        self.lattice_system = lattice_system

    def __repr__(self):
        """repr"""
        return f"SpaceGroup({self.number}, {self.symbol})"


def _load_spg_group_data():
    """ Does the loading, kept in function so temporary data is disposed"""
    # This will be used to look up the Bravais lattice definition

    spacegroup_number_to_lattice_system = \
        ["a" for _ in range(1, 3)] + \
        ["m" for _ in range(3, 16)] + \
        ["o" for _ in range(16, 75)] + \
        ["t" for _ in range(75, 143)] + \
        ["h" for _ in range(143, 195)] + \
        ["c" for _ in range(195, 231)]

    spacegroup_symbol_to_bravais_symbol = {
        "A": "S",
        "B": "S",
        "C": "S",
        "F": "F",
        "I": "I",
        "P": "P",
        "R": "R",
    }

    spacegroup_names = {}
    spacegroup_short_names = {}
    for i in range(1, 531):
        group_data = spglib.get_spacegroup_type(i)

        number = group_data.number

        spacegroup_names[number] = group_data.international_full
        spacegroup_short_names[number] = group_data.international_short

    # Make a lookup for spacegroups
    spacegroup_to_magnetic_group = defaultdict(list[int])
    for i in range(1,1652):
        metadata = spglib.get_magnetic_spacegroup_type(i)
        spacegroup_to_magnetic_group[metadata.number].append(i)

    # Create magnetic groups
    magnetic_groups = []
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
        magnetic_groups.append(group)



    # Create spacegroups
    lattice_symbol_to_spacegroups = defaultdict(list[SpaceGroup])
    spacegroups = []


    for i in range(1, 231):

        # Get the relevant data

        op_data = spglib.get_symmetry_from_database(i)
        corresponding_magnetic_groups = [magnetic_groups[idx-1] for idx in spacegroup_to_magnetic_group[i]]

        # Classify
        name = spacegroup_names[i]
        short_name = spacegroup_short_names[i]
        lattice_type_from_name = name[0]

        first_letter = spacegroup_number_to_lattice_system[i-1]
        second_letter = spacegroup_symbol_to_bravais_symbol[lattice_type_from_name]
        bravais_lattice = first_letter + second_letter

        # Load the operations
        translations = op_data["translations"]
        rotations = op_data["rotations"]

        operations = []
        for translation, rotation, in zip(translations, rotations):
            op = SpaceOperation.from_numpy(rotation, translation)
            operations.append(op)

        group = SpaceGroup(
            number=i,
            international_symbol=name,
            short_symbol=short_name,
            operations=operations,
            magnetic_variants=corresponding_magnetic_groups,
            lattice_system=lattice_system_letter_lookup[first_letter])

        spacegroups.append(group)
        lattice_symbol_to_spacegroups[bravais_lattice].append(group)

    return spacegroups, lattice_symbol_to_spacegroups, magnetic_groups

# Load the data
spacegroups, spacegroup_lattice_symbol_lookup, magnetic_groups = _load_spg_group_data()

# Lookup for spacegroups by long/short international symbol
spacegroup_symbol_lookup: dict[str, SpaceGroup] = {group.symbol: group for group in spacegroups}
spacegroup_symbol_lookup.update({group.short_symbol: group for group in spacegroups})

magnetic_group_symbol_lookup: dict[str, MagneticSpaceGroup] = {group.symbol: group for group in magnetic_groups}

# reverse the spacegroup/lattice relationship ????
spacegroup_symbol_to_lattice = {
    spacegroup.symbol: lattice_symbol
        for lattice_symbol in spacegroup_lattice_symbol_lookup
        for spacegroup in spacegroup_lattice_symbol_lookup[lattice_symbol]}

_lowercase_space_group_lookup: dict[str, SpaceGroup] = {
    key.lower(): spacegroup_symbol_lookup[key]
        for key in spacegroup_symbol_lookup}

_lowercase_magnetic_group_lookup: dict[str, MagneticSpaceGroup] = {
    key.lower(): magnetic_group_symbol_lookup[key]
        for key in magnetic_group_symbol_lookup}


class NoSuchGroup(Exception):
    """ Raised when a group is not found. """


def spacegroup_by_name(name: str) -> SpaceGroup:
    """ Get a spacegroup by name"""

    lower_name = name.lower()

    if lower_name in _lowercase_space_group_lookup:
        return _lowercase_space_group_lookup[lower_name]

    # Failure branch

    potential_names = {key.lower(): key for key in spacegroup_symbol_lookup.keys()}
    similar = get_close_matches(name, potential_names.keys(), n=3, cutoff=0.4)

    if similar:
        suggestion_string = ", ".join([f"'{potential_names[s]}'" for s in similar])
        message_string = f"Unknown space group '{name}', perhaps you meant {suggestion_string} or something similar."

    else:
        message_string = f"Unknown space group '{name}', doesn't seem to be even close to a magnetic space group."

    raise NoSuchGroup(message_string)


def magnetic_spacegroup_by_name(name: str) -> MagneticSpaceGroup:
    """ Get a magnetic spacegroup by name

    TODO: make this work more like the normal spacegroup one
    """

    lower_name = name.lower()
    if lower_name in _lowercase_magnetic_group_lookup:
        return _lowercase_magnetic_group_lookup[lower_name]

    similar = get_close_matches(lower_name, _lowercase_magnetic_group_lookup.keys(), n=3)
    formatted = [_lowercase_magnetic_group_lookup[key].symbol for key in similar]

    if similar:
        suggestion_string = ", ".join([f"'{s}'" for s in formatted])
        message_string = f"Unknown space group '{name}', perhaps you meant {suggestion_string} or something similar."

    else:

        message_string = f"Unknown space group '{name}', doesn't seem to be even close to a magnetic space group."

    raise NoSuchGroup(message_string)


if __name__ == "__main__":
    print(magnetic_group_symbol_lookup)
