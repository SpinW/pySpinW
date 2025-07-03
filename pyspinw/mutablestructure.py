""" TODO: WIP """

from pyspinw.couplinggroup import CouplingGroup
from pyspinw.gui.symmetry_settings import SymmetrySettings
from pyspinw.site import LatticeSite


class MutableStructure:
    """ A container for working magnetic structures

    A structure like this is important for synchronising the list of sites with the
    """

    def __init__(self,
                 sites: list[LatticeSite],
                 coupling_groups: list[CouplingGroup],
                 symmetry: SymmetrySettings):

        self._sites = sites
        self._coupling_groups = coupling_groups
        self.symmetry = symmetry


        self._implied_sites = []
        self._implied_site_to_site = []
        self._site_to_implied_site = []

        self._update_symmmetry_sites()

    @property
    def sites(self):
        """ Getter for sites """

    def _update_symmmetry_sites(self):
        """ Update the sites

        This will recalculate the implied sites and the arrays that say how they are referenced to each other
        """
        implied_sites = []
        implied_site_to_site = []
        site_to_implied_site = []

        count = 0
        for site_index, site in enumerate(self._sites):
            extra_sites = self.symmetry.magnetic_group.duplicates(site)

            end_count = count + len(extra_sites)

            implied_sites += extra_sites
            implied_site_to_site += [site_index for _ in extra_sites]

            site_to_implied_site.append([i for i in range(count, end_count)])

            count = end_count

        self._implied_sites = implied_sites
        self._implied_site_to_site = implied_site_to_site
        self._site_to_implied_site = site_to_implied_site

    def check_symmetry(self):
        """ Check the symmetry of this system and the system definition are consitent """

    def propose_symmetry(self):
        """ Propose symmetries that would be a good choice for this system """

    def _make_real(self, implied_site_index):
        """ Take an implied site and turn it into a real site"""

    def remove_site(self):
        """ Remove a site from the definition"""

        # Work out what couplings need to change / disappear

        # Remove site and update



    def remove_coupling_group(self):
        """ Remove a coupling group from the system definition """


