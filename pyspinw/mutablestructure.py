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

        self.sites = sites
        self.coupling_groups

    def update_symmmetry_sites(self):
        pass

    def remove_site(self):


        # Work out what couplings need to change / disappear

        # Remove site and update

        pass


    def remove_coupling_group(self):
        pass


    