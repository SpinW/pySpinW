import numpy as np

import spglib

from pyspinw.serialisation import SPWSerialisable
from pyspinw.site import ImpliedLatticeSite, LatticeSite
from pyspinw.symmetry.group import SymmetryGroup, database, ExactMatch
from pyspinw.symmetry.operations import MagneticOperation
from pyspinw.tolerances import tolerances


class MagneticSpaceGroup(SymmetryGroup, SPWSerialisable):
    """ Representation of a magnetic space group"""

    serialisation_name = "magnetic_group"

    def __init__(self, operations: list[MagneticOperation]):
        self.operations = operations

    def __repr__(self):
        """repr"""
        return f"MagneticSpaceGroup({self.number}, {self.symbol})"

    def implied_sites_for(self, site: LatticeSite) -> list[ImpliedLatticeSite]:
        """ Find "duplicate" sites of a given site """
        coordinates = site.values.reshape(1, -1) % 1

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

    def get_spacegroup(self):

        space_operations = [operation.space_operation().text_form for operation in self.operations]

        match_data = database.spacegroups_with_operations("; ".join(space_operations))

        if isinstance(match_data, ExactMatch):
            return match_data.spacegroup

        else:
            raise ValueError("Failed to find spacegroup exactly matching this magnetic group")

    @staticmethod
    def closure(operations: list[MagneticOperation], centerings: list[MagneticOperation]):
        # Find the closure of the group
        all_operations = set(centering.and_then(op)
                             for centering in centerings
                              for op in operations)
    
        last_operations_count = 0
        while last_operations_count < len(all_operations):
            last_operations_count = len(all_operations)

            new_operations = set(
                                op1.and_then(op2)
                                for op1 in all_operations
                                for op2 in all_operations)

            all_operations = all_operations.union(new_operations)

        return MagneticSpaceGroup(list(all_operations))


def load_magnetic_spacegroups():
    for uni in range(1, 1652):
        symmetry_data = spglib.get_magnetic_symmetry_from_database(uni)
        metadata = spglib.get_magnetic_spacegroup_type(uni)

        print(metadata)

if __name__ == "__main__":
    load_magnetic_spacegroups()