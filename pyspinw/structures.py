""" Magnetic structures """


import numpy as np

from pyspinw.serialisation import SPWSerialisable
from pyspinw.site import LatticeSite
from pyspinw.symmetry.group import SpaceGroup, MagneticSpaceGroup, SymmetryGroup, database
from pyspinw.symmetry.supercell import Supercell, TrivialSupercell
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.tolerances import tolerances
from pyspinw.util import connected_components, arraylike_equality


class Structure(SPWSerialisable):
    """ Representation of the magnetic structure """

    def __init__(self,
                 sites: list[LatticeSite],
                 unit_cell: UnitCell,
                 spacegroup: SymmetryGroup | None = None,
                 supercell: Supercell | None = None):

        self._input_sites = sites
        self._unit_cell = unit_cell


        self._spacegroup = database.spacegroups[1] if spacegroup is None else spacegroup
        self._supercell = TrivialSupercell() if supercell is None else supercell

        self._sites: list[LatticeSite] = self._extended_sites()

    def full_structure_site_list(self):
        """ All the sites in the structure"""
        all_sites = []

        # unit cell stuff
        cell_sites = self._extended_sites()

        # supercell
        for cell in self.supercell.cells():
            all_sites += [LatticeSite(
                i=site.i + cell.i,
                j=site.j + cell.j,
                k=site.k + cell.k,
                supercell_moments=site.moment_data,
                name=site.name,
                S=site.magnitude,
            ) for site in cell_sites]

        return all_sites

    def matplotlib_site_data(self):
        """ Data for making matplotlib scatter plots of sites """
        array_data = np.array([self.unit_cell.fractional_to_cartesian(site._ijk)
                         for site in self.full_structure_site_list()])

        return array_data[:,0], array_data[:, 1], array_data[:, 2]

    def _extended_sites(self)  -> list[LatticeSite]:
        """ All the sites, including those implied by symmetry """
        site_list = self._input_sites.copy()
        for site in self._input_sites:
            site_list += self._spacegroup.implied_sites_for(site)

        # Check for collisions, if there is an input site that
        # collides with an implied site, choose the input site
        # if there is a collision between two implied sites,
        # create a new one with:
        #  1: the same momement if the moments agree
        #  2: zero moment (of same shape) if they do not agree

        n_sites_raw = len(site_list)
        collisions = np.zeros((n_sites_raw, n_sites_raw), dtype=bool)
        for site_index_1, site_1 in enumerate(site_list):
            for site_index_2, site_2 in enumerate(site_list):
                # Do they have the same position
                collisions[site_index_1, site_index_2] = \
                    np.all(np.abs(site_1.ijk - site_2.ijk) < tolerances.SAME_SITE_ABS_TOL)

        # Get all the groups of collisions
        components = connected_components(collisions)

        input_site_ids = [site._unique_id for site in self._input_sites]

        # We want one element of each component, based on the criteria above
        unique_sites = []
        for component in components:
            sites = [site_list[i] for i in component]

            # Do we have a non-implied site in the list
            potential_choices = []
            for site in sites:
                if site._unique_id in input_site_ids:
                    potential_choices.append(site)

            if len(potential_choices) > 0:
                if len(potential_choices) == 1:
                    unique_sites.append(potential_choices[0])
                else:
                    raise ValueError("Multiple input sites at same spatial location")

                continue

            # At this point, there are no non-implied sites

            # Do they all have the same moment
            same_moment = True
            site_1 = sites[0]
            for site_2 in sites[1:]:
                if not arraylike_equality(
                          site_1._moment_data,
                          site_2._moment_data,
                          tolerances.SAME_SITE_ABS_TOL):

                    same_moment = False

            if same_moment:
                unique_sites.append(LatticeSite(
                    site_1.i,
                    site_1.j,
                    site_1.k,
                    supercell_moments=site_1.moment_data,
                    name = site_1.name, # TODO: Check if this is sensible
                    S=site_1.magnitude,
                ))

            else:
                unique_sites.append(LatticeSite(
                    site_1.i,
                    site_1.j,
                    site_1.k,
                    supercell_moments=np.zeros_like(site_1.moment_data),
                    name=site_1.name,  # TODO: Check if this is sensible
                    S=site_1.magnitude,
                ))

        return unique_sites

    def expansion_site_mapping(self):
        """ Expand supercell into a single, bigger cell """
        # Calculate new cell
        scale = self.supercell.cell_size()

        big_cell = self.unit_cell.updated(
            a=self.unit_cell.a * scale[0],
            b=self.unit_cell.b * scale[1],
            c=self.unit_cell.c * scale[2])

        # Create a mapping between sites and offsets to the new sites
        mapping: dict[tuple[LatticeSite, tuple[int, int, int]], LatticeSite] = {}
        for offset in self.supercell.cells():
            for site in self.sites:
                position = self.supercell.fractional_in_supercell(site.ijk, offset)
                moment = self.supercell.moment(site, cell_offset=offset)

                new_site = LatticeSite(
                    i=position[0],
                    j=position[1],
                    k=position[2],
                    supercell_moments=moment,
                    g=site.g,
                    name=site.name,
                    S=site.magnitude)

                mapping[(site._unique_id, offset.as_tuple)] = new_site

        return big_cell, mapping

    def expand(self):
        """ Expand supercell into a single, bigger cell """
        cell, mapping = self.expansion_site_mapping()


        return Structure(
            sites=[site for site in mapping.values()],
            unit_cell=cell,
            spacegroup=self.spacegroup.for_supercell(self.supercell),
            supercell=TrivialSupercell(scaling=(1,1,1))
        )



    def _build_sites(self):
        """ Updates the site variable on parameter changes"""
        self._sites = self._extended_sites()

    @property
    def sites(self) -> list[LatticeSite]:
        """ Get the sites used to define the structure (but implied by symmetry) """
        return self._input_sites.copy()

    @sites.setter
    def sites(self, sites):
        """ Set the input sites """
        self._input_sites = sites
        self._sites = self._build_sites()

    @property
    def spacegroup(self) -> SpaceGroup | MagneticSpaceGroup:
        """ Get the spacegroup"""
        return self._spacegroup

    @spacegroup.setter
    def spacegroup(self, spacegroup: SpaceGroup | MagneticSpaceGroup):
        """ Set the spacegroup"""
        self._spacegroup = spacegroup
        self._build_sites()

    @property
    def unit_cell(self) -> UnitCell:
        """ Get the unit cell"""
        return self._unit_cell

    @unit_cell.setter
    def unit_cell(self, unit_cell: UnitCell):
        """ Set the unit cell"""
        self._unit_cell = unit_cell
        self._build_sites()

    @property
    def supercell(self) -> Supercell:
        """ Get the supercell """
        return self._supercell

    @supercell.setter
    def supercell(self, supercell: Supercell):
        """ Set the supercell """
        self._supercell = supercell
        self._build_sites()
