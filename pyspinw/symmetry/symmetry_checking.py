""" Checking of consistency symmetry """
from collections import defaultdict

import numpy as np

from pyspinw.site import LatticeSite
from pyspinw.symmetry.group import MagneticSpaceGroup
from pyspinw.symmetry.supercell import Supercell
from pyspinw.tolerances import tolerances


def check_individual_moment_consistency(position: np.ndarray,
                                        moment: np.ndarray,
                                        group: MagneticSpaceGroup,
                                        verbose: bool=False) -> list[str]:

    position_and_moment = np.concatenate((position, moment)).reshape(-1, 6)
    info = []
    for symmetry in group.operations:
        transformed = symmetry(position_and_moment)

        # If the location is not invariant with respect to the symmetry operation, then it doesn't matter, continue
        new_position = transformed[0, :3]
        if np.any(np.abs(new_position - position) > tolerances.SAME_SITE_ABS_TOL):
            if verbose:
                info.append(f"Site at ({position[0]}, {position[1]}, {position[2]}) is NOT "
                            f"INVARIANT under '{symmetry.text_form}'")

            continue

        # If the location is invariant, the moment should be too
        new_moment = transformed[0, 3:]
        if np.any(np.abs(new_moment - moment) > tolerances.SAME_SITE_ABS_TOL):

            info.append(f"Site at ({position[0]}, {position[1]}, {position[2]}) "
                            f"with moment ({moment[0]}, {moment[1]}, {moment[2]}) "
                            f"is inconsistent with symmetry operation '{symmetry.text_form}'")

        else:
            if verbose:
                info.append(f"Site at ({position[0]}, {position[1]}, {position[2]}) "
                            f"with moment ({moment[0]}, {moment[1]}, {moment[2]}) "
                            f"is CONSISTENT with symmetry operation '{symmetry.text_form}'")

    return info


def check_supercell_moment_consistency(
        supercell: Supercell,
        group: MagneticSpaceGroup,
        sites: list[LatticeSite],
        verbose: bool=False):
    """ Check consistency of magnetic moments with the spacegroup"""
    info = defaultdict(list[str])

    for cell_offset in supercell.cells():
        for site in sites:
            moment = supercell.moment(site, cell_offset)
            cell_site_errors = check_individual_moment_consistency(site.ijk, moment, group, verbose=verbose)
            if cell_site_errors:
                info[cell_offset.as_tuple] += cell_site_errors

    return info