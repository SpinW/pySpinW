""" Checking of consistency symmetry """
from collections import defaultdict

import numpy as np

from pyspinw.coupling import Coupling
from pyspinw.site import LatticeSite
from pyspinw.symmetry.group import MagneticSpaceGroup
from pyspinw.symmetry.supercell import Supercell
from pyspinw.tolerances import tolerances


def check_supercell_moment_consistency(
        supercell: Supercell,
        group: MagneticSpaceGroup,
        sites: list[LatticeSite]):
    """ Check consistency of magnetic moments in a supercell with the magnetic spacegroup"""
    position_and_moments = []
    offsets = []
    for cell_offset in supercell.cells():
        for site in sites:
            position_and_moments.append(supercell.cell_position_and_moment(site, cell_offset))
            offsets.append(cell_offset.as_tuple)

    position_and_moments = np.array(position_and_moments)
    offsets = np.array(offsets)

    info = defaultdict(list[str])
    for symmetry in group.operations:
        transformed = symmetry(position_and_moments)

        # We want to search for position invariants that do not leave the momentum unchanged
        diffs = np.abs(position_and_moments - transformed)
        same_positions = np.all(diffs[:, :3] < tolerances.SAME_SITE_ABS_TOL, axis=1)
        same_momenta = np.all(diffs[:, 3:] < tolerances.SAME_SITE_ABS_TOL, axis=1)

        problems = same_positions & ~same_momenta
        problem_sites = position_and_moments[problems, :]
        problem_offsets = offsets[problems, :]

        for site, raw_offset in zip(problem_sites, problem_offsets):
            offset = tuple(int(x) for x in raw_offset)
            info[offset].append(f"Site at ({site[0]:.4g}, {site[1]:.4g}, {site[2]:.4g}) "
                                f"with moment ({site[3]:.4f}, {site[4]:.4f}, {site[5]:.4f}) "
                                f"in cell at ({offset[0]}, {offset[1]}, {offset[2]}) "
                                f"is magnetically inconsistent under '{symmetry.text_form}'")

    return info

def check_coupling_consistency(sites: list[LatticeSite], couplings: list[Coupling]):
    """ Check that a coupling actually does something, not cancelled by symmetry """
    for coupling in couplings:

        # Are they referring to the same site
        if coupling.site_1.parent_site._unique_id == coupling.site_2.parent_site._unique_id:
            # We want to check whether (R1 S)^T M (R2 S) is constant
            #  As the magnitude of S can be different when actually running the calculation,
            #  this can only happen when the constant is zero,
            #  which is when R1 M R2 is antisymmetric

            pass


