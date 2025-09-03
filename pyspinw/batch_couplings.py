""" Functions and classes for automatically creating multiple couplings"""

from collections import defaultdict

import re
from dataclasses import dataclass

import numpy as np

from pyspinw.lattice_distances import find_relative_positions
from pyspinw.cell_offsets import CellOffset
from pyspinw.site import LatticeSite, ImpliedLatticeSite
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.tolerances import tolerances


@dataclass
class AbstractCoupling:
    """ Used to internally represent some coupling properties, and some properties that are only needed internally"""

    name: str
    site_1: LatticeSite
    site_2: LatticeSite
    cell_offset: CellOffset
    distance: float
    order: int

def approx_equal(v1, v2):
    """ Are two vectors approximately the same (within tolerance)"""
    return np.all(np.abs(v1 - v2) < tolerances.VECTOR_TOL)

def approx_equal_direction(v1, v2):
    "Are two vectors in the same direction (within tolerance)"
    return np.all(np.abs(v1 - v2) < tolerances.VECTOR_TOL) or np.all(np.abs(v1 + v2) < tolerances.VECTOR_TOL)


def apply_naming_convention(naming_pattern: str,
                            name_1: str,
                            name_2: str,
                            type_symbol: str,
                            order: int,
                            xyz_direction: np.ndarray):
    """Convert a naming string into an actual name

    :param naming_pattern: string describing how to name things
    :param name_1: name of the first site - specified with $SITE1$
    :param name_2: name of the second site - specified with $SITE2$
    :param type_symbol: string that denotes the type of the coupling (e.g. 'J', 'DM')
    :param order: index of the "shell", i.e. an index that increases with distance - use $ORDER$
    :param xyz_direction: direction of coupling, use $DIRECTION$ to give a string that denotes
                           this concisely (hopefully)

    """
    if "[direction]" in naming_pattern:

        if np.all(np.abs(xyz_direction) < 1e-10):
            direction_string = ""

        else:
            normalised = xyz_direction / np.max(np.abs(xyz_direction))

            if approx_equal_direction(normalised, [1,0,0]):
                direction_string = "x"
            elif approx_equal_direction(normalised, [0,1,0]):
                direction_string = "y"
            elif approx_equal_direction(normalised, [0,0,1]):
                direction_string = "z"
            else:
                direction_string = f"[{normalised[0]:.3g},{normalised[1]:.3g},{normalised[2]:.3g}]"

        naming_pattern = re.sub(r"\[direction]", direction_string, naming_pattern)

    # Less complicated substitutions
    naming_pattern = re.sub(r"\[site1]", name_1, naming_pattern)
    naming_pattern = re.sub(r"\[site2]", name_2, naming_pattern)
    naming_pattern = re.sub(r"\[order]", str(order), naming_pattern)
    naming_pattern = re.sub(r"\[type]", type_symbol, naming_pattern)

    return naming_pattern



default_naming_pattern = "[type][order]:[site1]-[site2]" #[type]_[order]([site1], [site2])"

def batch_couplings(sites: list[LatticeSite],
                    unit_cell: UnitCell,
                    max_distance: float,
                    naming_pattern: str=default_naming_pattern,
                    type_symbol: str="J",
                    both_directions: bool=False):
    """ Find all the couplings within a certain distance

    :param sites: List of LatticeSite or ImpliedLatticeSite to find couplings between
    :param unit_cell: Unit cell, needed for working out cartesian distances
    :param max_distance: Maximum distance to get couplings for
    :param naming_convention: Formatting string for the naming, see `apply_naming_convention` for details
    :param both_directions: Couplings in both directions


    Most of the work here is about constraining by order and giving it a sensible name
    """
    unique_sites = set([site.parent_site if isinstance(site, ImpliedLatticeSite) else site for site in sites])
    unique_index = {site: index for index, site in enumerate(unique_sites)}

    # Get all the distances we'll need

    pair_data = defaultdict(dict)
    for site_1_index, site_1 in enumerate(sites):

        site_2_start_index = 0 if both_directions else site_1_index

        for site_2 in sites[site_2_start_index:]:

            allow_self = site_1 is not site_2 # Flag to include (0,0,0) offset

            relative_position = site_2.ijk - site_1.ijk
            same_cell_relative_position = (site_2.ijk - site_1.ijk) % 1.0
            cell_correction = np.array(relative_position - same_cell_relative_position, dtype=float)


            positions = find_relative_positions(
                same_cell_relative_position,
                unit_cell._xyz,
                max_distance=max_distance,
                allow_self=allow_self)

            pair_data[site_1][site_2] = [(site_1, site_2, vector, cell_offset - cell_correction, distance)
                                         for cell_offset, vector, distance in positions.expand()]


    # Reduce down to unique sites for the purpose of naming and for finding shells
    reduced_pair_data = defaultdict(lambda: defaultdict(list))
    for site_1 in pair_data:
        for site_2 in pair_data[site_1]:

            if isinstance(site_1, ImpliedLatticeSite):
                root_site_1 = site_1.parent_site
            else:
                root_site_1 = site_1

            if isinstance(site_2, ImpliedLatticeSite):
                root_site_2 = site_2.parent_site
            else:
                root_site_2 = site_2


            # Sort so that we have a unique ordering of the pairs
            if unique_index[root_site_1] < unique_index[root_site_2]:
                pair = root_site_1, root_site_2
            else:
                pair = root_site_2, root_site_1

            reduced_pair_data[pair[0]][pair[1]] += pair_data[site_1][site_2]


    all_couplings = []

    # Find shells and create abstract couplings
    for root_site_1 in reduced_pair_data:
        for root_site_2 in reduced_pair_data[root_site_1]:
            all_links = reduced_pair_data[root_site_1][root_site_2]

            if not all_links:
                continue

            distances = sorted([(i, data[4]) for i, data in enumerate(all_links)], key=lambda x: x[1])

            # Assign indices to distances

            orders = [(distances[0][0], 1)] # Pairs of (index, order)

            if len(distances) > 1:
                last_distance = distances[0][1]
                order = 1
                for index, distance in distances[1:]:
                    if distance - last_distance > tolerances.COUPLING_ORDER_THRESHOLD:
                        order += 1

                    orders.append((index, order))
                    last_distance = distance

            # Create object
            for link_index, order in orders:
                site_1, site_2, vector, cell_offset_raw, distance = all_links[link_index]

                if isinstance(site_1, ImpliedLatticeSite):
                    site_1_name = site_1.parent_site.name
                else:
                    site_1_name = site_1.name

                if isinstance(site_2, ImpliedLatticeSite):
                    site_2_name = site_2.parent_site.name
                else:
                    site_2_name = site_2.name


                cell_offset = CellOffset(i=cell_offset_raw[0], j=cell_offset_raw[1], k=cell_offset_raw[2])

                name = apply_naming_convention(naming_pattern,
                                               site_1_name,
                                               site_2_name,
                                               type_symbol,
                                               order,
                                               vector)

                all_couplings.append(AbstractCoupling(
                    name = name,
                    site_1 = site_1,
                    site_2 = site_2,
                    cell_offset = cell_offset,
                    distance = distance,
                    order = order ))

    # Return sorted list by distance
    return sorted(all_couplings, key=lambda coupling: coupling.distance)


