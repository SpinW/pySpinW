""" TODO: WIP """
from dataclasses import dataclass
from enum import Enum

import numpy as np

from pyspinw.coupling import Coupling
from pyspinw.couplinggroup import CouplingGroup
from pyspinw.gui.symmetry_settings import SymmetrySettings, DEFAULT_SYMMETRY
from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.tolerances import tolerances


@dataclass
class BoundCouplingGroup:
    """ This links a coupling group to sites"""

    def __init__(self, coupling_group: CouplingGroup, indices: list[int]):
        self.coupling_group = coupling_group
        self.indices = indices

    def couplings(self, sites: list[LatticeSite], unit_cell: UnitCell) -> list[Coupling]:
        """ Return couplings for the coupling group """
        return self.coupling_group.couplings([sites[index] for index in self.indices], unit_cell=unit_cell)

class SymmetryConflictType(Enum):
    """ Types of symmetry conflict"""

    POSITION_CONFLICT = "position conflict"
    MAGNETIC_CONFLICT = "magnetic conflict"
    COUPLING_CANCELLED = "coupling cancelled"

@dataclass
class SymmetryConflict:
    """ Representation of a symmetry conflict"""

    conflict_type: SymmetryConflictType
    message: str

class MutableStructure:
    """ A container for working magnetic structures

    A structure like this is important for synchronising the list of sites with the
    """

    def __init__(self,
                 sites: list[LatticeSite] | None = None,
                 coupling_groups: list[CouplingGroup] | None = None,
                 symmetry: SymmetrySettings = DEFAULT_SYMMETRY):

        self._sites = [] if sites is None else sites
        self._coupling_groups = [] if coupling_groups is None else coupling_groups
        self._symmetry = symmetry

        self._implied_sites = []
        self._implied_site_to_site = []
        self._site_to_implied_site = []

        self._update_symmmetry_sites()

    @property
    def sites(self):
        """ Getter for sites """
        return self._sites

    def _update_symmmetry_sites(self):
        """ Update the sites

        This will recalculate the implied sites and the arrays that say how they are referenced to each other
        """
        implied_sites = []
        implied_site_to_site = []
        site_to_implied_site = []

        count = 0
        for site_index, site in enumerate(self._sites):
            extra_sites = self.symmetry.magnetic_group.implied_sites_for(site)

            end_count = count + len(extra_sites)

            implied_sites += extra_sites
            implied_site_to_site += [site_index for _ in extra_sites]

            site_to_implied_site.append([i for i in range(count, end_count)])

            count = end_count

        self._implied_sites = implied_sites
        self._implied_site_to_site = implied_site_to_site
        self._site_to_implied_site = site_to_implied_site

    def _all_sites(self):
        return self._sites + self._implied_sites

    def _couplings(self):

        couplings = []
        all_sites = self._all_sites()

        for coupling_group in self._coupling_groups:
            couplings += coupling_group.couplings(all_sites, self._symmetry.unit_cell)

    def check_symmetry(self) -> list[SymmetryConflict]:
        """ Check the symmetry of this system and the system definition are consistent """
        check_results = []

        # Check site symmetry
        all_sites = self._all_sites()
        for i, site1 in enumerate(all_sites):
            for site2 in all_sites[:i]:
                if np.all(np.abs(site1.ijk - site2.ijk) < tolerances.SAME_SITE_ABS_TOL):

                    if np.all(np.abs(site1.m - site2.m) > tolerances.SAME_SITE_ABS_TOL):
                        conflict = SymmetryConflict(
                            SymmetryConflictType.POSITION_CONFLICT,
                                f"{site1} and {site2} share the same position and moment")

                    else:
                        conflict = SymmetryConflict(
                            SymmetryConflictType.MAGNETIC_CONFLICT,
                                f"{site1} and {site2} share the same position, but with different moments")

                    check_results.append(conflict)

        # Check coupling symmetry
        couplings = self._couplings()

        # TODO

        return check_results

    def propose_symmetry(self):
        """ Propose symmetries that would be a good choice for this system """

    def _make_real(self, implied_site_indices: list[int]):
        """ Take an implied site and turn it into a real site"""
        new_sites = [self._implied_sites[i].reify() for i in implied_site_indices]

        self._sites += new_sites

        self._update_symmmetry_sites()




    def remove_site(self):
        """ Remove a site from the definition"""

        # Work out what couplings need to change / disappear

        # Remove site and update



    def remove_coupling_group(self):
        """ Remove a coupling group from the system definition """


