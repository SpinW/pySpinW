from collections import defaultdict

import numpy as np
import spglib

from pyspinw.site import LatticeSite
from pyspinw.symmetry.operations import MagneticOperation, SpaceOperation
from pyspinw.symmetry.data.msg_symbols import msg_symbols
from pyspinw.tolerances import tolerances


class SymmetryGroup:
    pass

class MagneticSpaceGroup(SymmetryGroup):
    def __init__(self, number: int, symbol: str, operations: list[MagneticOperation]):
        self.number = number
        self.symbol = symbol
        self.operations = operations

    def __repr__(self):
        return f"SpaceGroup({self.number}, {self.symbol})"

    def duplicates(self, site: LatticeSite):
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
            new_site = LatticeSite.from_coordinates(ijkm.reshape(-1), name=site.name + f" [{i+1}]")

            new_sites.append(new_site)

        return new_sites



class SpaceGroup(SymmetryGroup):
    def __init__(self, number, international_symbol, operations, magnetic_variants: list[MagneticSpaceGroup]):
        self.number = number
        self.symbol = international_symbol
        self.operations = operations
        self.magnetic_variants = magnetic_variants

    def __repr__(self):
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
    for i in range(1, 531):
        group_data = spglib.get_spacegroup_type(i)
        number = group_data["number"]
        name = group_data["international_full"]
        spacegroup_names[number] = name

    # Make a lookup for spacegroups
    spacegroup_to_magnetic_group = defaultdict(list[int])
    for i in range(1,1652):
        metadata = spglib.get_magnetic_spacegroup_type(i)
        spacegroup_to_magnetic_group[metadata["number"]].append(i)

    # Create magnetic groups
    magnetic_groups = []
    for i in range(1,1652):

        op_data = spglib.get_magnetic_symmetry_from_database(i)
        metadata = spglib.get_magnetic_spacegroup_type(i)
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
            international_symbol = name,
            operations = operations,
            magnetic_variants=corresponding_magnetic_groups)

        spacegroups.append(group)
        lattice_symbol_to_spacegroups[bravais_lattice].append(group)

    return spacegroups, lattice_symbol_to_spacegroups, magnetic_groups

# Load the data
spacegroups, spacegroup_lattice_symbol_lookup, magnetic_groups = _load_spg_group_data()
spacegroup_symbol_lookup = {group.symbol: group for group in spacegroups}
magnetic_group_symbol_lookup = {group.symbol: group for group in magnetic_groups}

if __name__ == "__main__":
    print(spacegroup_lattice_symbol_lookup)