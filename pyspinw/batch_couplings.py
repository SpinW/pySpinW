from collections import defaultdict

import re

import numpy as np

from lattice_distances import find_relative_positions
from pyspinw.gui.cell_offsets import CellOffset
from pyspinw.site import LatticeSite, ImpliedLatticeSite
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.tolerances import tolerances


@dataclass
class AbstractCoupling:
    name: str
    site_1: LatticeSite
    site_2: LatticeSite
    cell_offset: CellOffset

def approx_equal(v1, v2):
    return np.all(np.abs(v1 - v2) < tolerances.VECTOR_TOL)

def approx_equal_direction(v1, v2):
        return np.all(np.abs(v1 - v2) < tolerances.VECTOR_TOL) or np.all(np.abs(v1 + v2) < tolerances.VECTOR_TOL)


def apply_naming_convention(naming_pattern: str, name_1: str, name_2: str, order: int, xyz_direction: np.ndarray, ijk_direction: np.ndarray):
    """Convert a naming string into an actual name

    :param naming_pattern: string describing how to name things
    :param name_1: name of the first site - specified with $SITE1$
    :param name_2: name of the second site - specified with $SITE2$
    :param order: index of the "shell", i.e. an index that increases with distance - use $ORDER$
    :param direction: direction of coupling, use $DIRECTION$ to give a string that denotes this consisely (hopefully)
    :param ijk_direction: direction in lattice

    """

    if "$DIRECTION$" in naming_pattern:
        # Gives xyz labels in preference to ijk
        if approx_equal_direction(xyz_direction, [1,0,0]):
            direction_string = "x"
        elif approx_equal_direction(xyz_direction, [0,1,0]):
            direction_string = "y"
        elif approx_equal_direction(xyz_direction, [0,0,1]):
            direction_string = "z"
        elif approx_equal_direction(ijk_direction, [1,0,0]):
            direction_string = "h"
        elif approx_equal_direction(ijk_direction, [0,1,0]):
            direction_string = "k"
        elif approx_equal_direction(ijk_direction, [0,0,1]):
            direction_string = "l"
        else:
            normalised = ijk_direction / np.max(ijk_direction)
            direction_string = f"{normalised[0]:3f},{normalised[1]:3f},{normalised[2]:3f}"

        naming_pattern = re.sub(r"\$DIRECTION\$", direction_string, naming_pattern)

    # Less complicated substitutions
    naming_pattern = re.sub(r"\$SITE1\$", name_1, naming_pattern)
    naming_pattern = re.sub(r"\$SITE2\$", name_2, naming_pattern)
    naming_pattern = re.sub(r"\$ORDER\$", order, naming_pattern)

    return naming_pattern


def batch_couplings(sites: list[LatticeSite], unit_cell: UnitCell, max_distance: float, naming_convention: str):
    """ Find all the couplings within a certain distance

    :param sites: List of LatticeSite or ImpliedLatticeSite to find couplings between
    :param unit_cell: Unit cell, needed for working out cartesian distances
    :param max_distance: Maximum distance to get couplings for
    :param naming_convention: Formatting string for the naming, see `apply_naming_convention` for details


    Most of the work here is about giving it a sensible name
    """


    unique_sites = set([site.parent_site if isinstance(site, ImpliedLatticeSite) else site for site in sites])
    unique_index = {site: index for index, site in enumerate(unique_sites)}

    # Get all the distances we'll need

    pair_data = defaultdict(dict)
    for site_1 in sites:
        for site_2 in sites:

            allow_self = site_1 is not site_2

            relative_base_distance = (site_2.ijk - site_1.ijk) % 1.0


            positions = find_relative_positions(relative_base_distance, unit_cell._xyz, max_distance=max_distance, allow_self=allow_self)
            pair_data[site_1][site_2] = [(site_1, site_2, vector, cell_offset, distance)
                                         for vector, cell_offset, distance in positions.expand()]


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


    # Find shells and create abstract couplings
    for root_site_1 in reduced_pair_data:
        for root_site_2 in reduced_pair_data[root_site_1]:
            all_links = reduced_pair_data[root_site_1][root_site_2]

            distances = sorted([(i, data[4]) for i, data in enumerate(all_links)], key=lambda x: x[1])


