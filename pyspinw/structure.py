""" Magnetic structures """
import logging

import numpy as np
from ase.data import chemical_symbols
from numpy._typing import ArrayLike

from pyspinw.cell_offsets import CellOffset
from pyspinw.exchangegroup import DirectionalityFilter
from pyspinw.lattice_distances import full_search_space
from pyspinw.serialisation import SPWSerialisable
from pyspinw.site import LatticeSite, ImpliedLatticeSite
from pyspinw.symmetry.group import SpaceGroup, MagneticSpaceGroup, SymmetryGroup, database
from pyspinw.symmetry.operations import SpaceOperation
from pyspinw.symmetry.supercell import Supercell, TiledSupercell
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.tolerances import tolerances
from pyspinw.util import connected_components, arraylike_equality, IncrementalPointHistogram

logger = logging.getLogger("Structure")

class Structure(SPWSerialisable):
    """ Representation of the magnetic structure """

    def __init__(self,
                 sites: list[LatticeSite],
                 unit_cell: UnitCell,
                 spacegroup: SymmetryGroup | None = None,
                 supercell: Supercell | None = None):

        self._input_sites = sites
        self._input_uid_to_site = {site.unique_id: site for site in sites}
        self._unit_cell = unit_cell

        self._spacegroup = database.spacegroups[0] if spacegroup is None else spacegroup
        self._supercell = TiledSupercell() if supercell is None else supercell

        self._sites: list[LatticeSite] = self._extended_sites()

        # Check that supercell components match site dimensions
        bad_sites = []
        for site in self.sites:
            if site.n_components() != self._supercell.n_components():
                bad_sites.append(site)

        if bad_sites:
            raise ValueError("Expected the shape of site spin data to match what the supercell requires "
                             f"({supercell.n_components()}-by-3), "
                             "bad sites are: " + ", ".join([site.name for site in bad_sites]))

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
                supercell_spins=site.spin_data,
                name=site.name
            ) for site in cell_sites]

        return all_sites

    def matplotlib_site_data(self):
        """ Data for making matplotlib scatter plots of sites """
        array_data = np.array([self.unit_cell.lattice_units_to_cartesian(site._ijk)
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
        #  1: the same momement if the spins agree
        #  2: zero spin (of same shape) if they do not agree

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

            # Do they all have the same spin
            same_spin = True
            site_1 = sites[0]
            for site_2 in sites[1:]:
                if not arraylike_equality(
                          site_1._spin_data,
                          site_2._spin_data,
                          tolerances.SAME_SITE_ABS_TOL):

                    same_spin = False

            if same_spin:
                unique_sites.append(LatticeSite(
                    site_1.i,
                    site_1.j,
                    site_1.k,
                    supercell_spins=site_1.spin_data,
                    name = site_1.name,
                    metadata = site_1.metadata.copy()
                ))

            else:
                unique_sites.append(LatticeSite(
                    site_1.i,
                    site_1.j,
                    site_1.k,
                    supercell_spins=np.zeros_like(site_1.spin_data),
                    name=site_1.name,
                    metadata=site_1.metadata.copy()
                ))

        return unique_sites

    def _expansion_site_mapping(self):
        """ Expand supercell into a single, bigger cell """
        # Calculate new cell
        scale = self.supercell.cell_size()

        big_cell = self.unit_cell.updated(
            a=self.unit_cell.a * scale[0],
            b=self.unit_cell.b * scale[1],
            c=self.unit_cell.c * scale[2])

        # Create a mapping between sites and offsets to the new sites
        mapping: dict[tuple[int, tuple[int, int, int]], LatticeSite] = {}
        for index, offset in enumerate(self.supercell.cells()):
            for site in self.sites:
                position = self.supercell.fractional_in_supercell(site.ijk, offset)
                spin = self.supercell.spin(site, cell_offset=offset)

                new_site = LatticeSite(
                    i=position[0],
                    j=position[1],
                    k=position[2],
                    supercell_spins=spin,
                    g=site.g,
                    name=f"{site.name}[{index}]",
                    metadata=site.metadata.copy())

                mapping[(site._unique_id, offset.as_tuple)] = new_site

        return big_cell, mapping

    def expand(self):
        """ Expand supercell into a single, bigger cell """
        cell, mapping = self._expansion_site_mapping()


        return Structure(
            sites=[site for site in mapping.values()],
            unit_cell=cell,
            spacegroup=self.spacegroup.for_supercell(self.supercell),
            supercell=TiledSupercell(scaling=(1, 1, 1))
        )

    def _without_nonmagnetic_with_removed_uids(self) -> tuple["Structure", set[int]]:
        """ Get a copy of the structure without magnetic sites, along with a list of uids for the removed sites """
        new_sites: list[LatticeSite] = []
        excluded_uids: set[int] = set()
        for site in self._input_sites:
            if np.allclose(site.spin_data, 0.0):
                excluded_uids.add(site.unique_id)
            else:
                new_sites.append(site)

        return Structure(new_sites, self.unit_cell, self.spacegroup, self.supercell), excluded_uids

    def without_nonmagnetic(self) -> "Structure":
        """ Get a copy of this structure with non-magnetic sites removed"""
        structure, _ = self._without_nonmagnetic_with_removed_uids()

        return structure

    def _build_sites(self):
        """ Updates the site variable on parameter changes"""
        self._sites = self._extended_sites()

    @property
    def sites(self) -> list[LatticeSite]:
        """ Get the sites used to define the structure and implied by symmetry """
        return self._sites.copy()
        # return self._input_sites.copy()

    def sites_by_name(self, name) -> list[LatticeSite]:
        """ Get sites where name matches regex"""
        return [site for site in self._sites if site.name.lower().startswith(name.lower())]

    def site_by_name(self, name: str):
        """ Get a single site by its name"""
        found = self.sites_by_name(name)
        if len(found) == 0:
            raise ValueError(f"No site matching '{name}' found")
        elif len(found) > 1:
            # Are any an exact match
            exact_matches = [site for site in found if site.name == name]

            if len(exact_matches) == 0:

                names = ", ".join([f"'{site.name}'" for site in found])

                raise ValueError(f"Multiple sites close to '{name}' found, but no exact match. "
                                 f"Close matches are {names}")

            elif len(exact_matches) == 1:
                return exact_matches[0]

            else:
                raise ValueError(f"Multiple sites exactly matching '{name}' found")

        else:
            return found[0]

    def sites_by_element(self, element: str | None) -> list[LatticeSite]:
        """ Get list of sites with the specified element in the metadata"""
        if element is not None and element not in chemical_symbols[1:]:
            raise ValueError(f"{element} is not a known element")
        return [site for site in self._sites if site.metadata.element == element]

    def without_elements(self, *elements):
        """ Get the structure without sites with given elements"""
        for element in elements:
            if element is not None and element not in chemical_symbols[1:]:
                raise ValueError(f"{element} is not a known element")

        new_sites = [site for site in self._sites if site.metadata.element not in elements]

        return Structure(new_sites, unit_cell=self.unit_cell, spacegroup=self.spacegroup, supercell=self.supercell)


    @property
    def text_summary(self) -> str:
        """ Textual details of this structure """
        lines = []
        lines.append(f"Unit Cell: {self.unit_cell.text_summary}")
        lines.append(f"Spacegroup: {self.spacegroup.preferred_symbol}")
        supercell_text_data = self.supercell.text_data()
        lines.append(f"Supercell: {supercell_text_data[0]}")
        lines += ["  " + s for s in supercell_text_data[1:]]
        lines.append("Sites:")
        for site in self._sites:

            is_not_input_chr = "" if site.unique_id in self._input_uid_to_site else "* "

            lines.append(f"  {is_not_input_chr}{site}")

        return "\n".join(lines)

    def print_summary(self):
        """ Print out details of this structure """
        print(self.text_summary)

    def all_neighbours(self,
                       site: LatticeSite,
                       n=1,
                       element: str | None = None,
                       direction_filter: DirectionalityFilter | None = None) \
            -> list[list[tuple[LatticeSite, CellOffset]]]:
        """ Get list-of-lists containing (i<=n)-th nearest neighbours along with their cell offsets

        :param n: maximum n for n-th nearest neighbour
        :param element: Restrict search to this element
        :param direction_filter: Filter the results using this

        :returns: a list of (site, offset) pairs for each order up-to and including `neighbour`
        """
        if element is not None and element not in chemical_symbols[1:]:
            raise ValueError(f"{element} is not an element")

        # Get a list of sites we want to check
        sites_to_check: list[LatticeSite] = []
        for test_site in self.sites:

            # Filter out sites without the specified element
            if element is not None and test_site.metadata.element != element:
                continue


            sites_to_check.append(test_site)

        # We want the n-th nearest neighbour, but we generate offsets in terms of distance
        #  so we need a strict bound on the distance that will cover all the points of a given order

        # possibly a very loose bound for systems with lots of sites, because it is based on
        #  having one site per unit cell
        search_radius = (n-1) * min([self.unit_cell.a, self.unit_cell.b, self.unit_cell.b])

        offsets = full_search_space(self.unit_cell._xyz, search_radius)

        n_checks = offsets.shape[0]*len(sites_to_check)

        if n_checks > 10_000:
            logger.warning(f"Neighbours is not optimised for large neighbour distances."
                           f"Finding the order-{n} neighbours entails checking {n_checks} "
                           f"unit cells")

        hist = IncrementalPointHistogram()

        for test_site in sites_to_check:

            lattice_vectors = test_site.ijk - site.ijk + offsets
            xyz_vectors = self.unit_cell.lattice_units_to_cartesian(lattice_vectors)
            distances = np.sqrt(np.sum(xyz_vectors**2, axis=1))

            for i in range(len(distances)):
                hist.add(distances[i], (test_site, xyz_vectors[i,:], offsets[i,:]))

        # Should have at least `neighbour+1` groups
        groups = hist.groups()[:n + 1]

        if direction_filter is None:
            return [[(test_site, CellOffset.coerce(offset))
                     for test_site, _, offset in order]
                        for order in groups]
        else:
            return [[(test_site, CellOffset.coerce(offset))
                     for test_site, vector, offset in order if direction_filter.accept(vector)]
                    for order in groups]



    def neighbours(self,
                   site: LatticeSite,
                   n=1,
                   element: str | None = None,
                   direction_filter: DirectionalityFilter | None = None):
        """ Get a list of n-th nearest neighbours, along with cell offsets

        :param n: maximum n for n-th nearest neighbour
        :param element: Restrict search to this element
        :param direction_filter: Filter the results using this

        :returns: a list of (site, offset) for n-th nearest neighbour
        """
        all_neighbors = self.all_neighbours(
                site = site,
                n=n,
                element=element,
                direction_filter=direction_filter)

        return all_neighbors[-1]

    def one_neighbour(self,
                      site: LatticeSite,
                      n=1,
                      element: str | None = None,
                      direction_filter: DirectionalityFilter | None = None):
        """ Get one, arbitrary, n-th nearest neighbours, along with cell offsets

        :param n: maximum n for n-th nearest neighbour
        :param element: Restrict search to this element
        :param direction_filter: Filter the results using this

        :returns: a list of (site, offset) for n-th nearest neighbour
        """
        neighbors = self.neighbours(
            site=site,
            n=n,
            element=element,
            direction_filter=direction_filter)

        if len(neighbors) == 0:
            raise ValueError("No neighbours that pass the filter")

        return neighbors[0]


    def symmetry_related(self, site: LatticeSite) -> list[tuple[LatticeSite, set[SpaceOperation]]]:
        """ Get a list of sites related to the specified site by symmetry, including the original site"""
        symmetry_related = []
        for other_site in self._sites:

            ops = self.spacegroup.operations_between_sites(site, other_site)
            if len(ops) > 0:
                symmetry_related.append((other_site, set(ops)))

        return symmetry_related

    @property
    def text_summary(self) -> str:
        """ Textual details of this structure """
        lines = []
        lines.append(f"Unit Cell: {self.unit_cell.text_summary}")
        lines.append(f"Spacegroup: {self.spacegroup.preferred_symbol}")
        supercell_text_data = self.supercell.text_data()
        lines.append(f"Supercell: {supercell_text_data[0]}")
        lines += ["  " + s for s in supercell_text_data[1:]]
        lines.append("Sites:")
        for site in self._sites:

            is_not_input_chr = "" if site.unique_id in self._input_uid_to_site else "* "

            lines.append(f"  {is_not_input_chr}{site}")

        return "\n".join(lines)

    def print_summary(self):
        """ Print out details of this structure """
        print(self.text_summary)

    @property
    def spacegroup(self) -> SpaceGroup | MagneticSpaceGroup:
        """ Get the spacegroup"""
        return self._spacegroup

    # @spacegroup.setter
    # def spacegroup(self, spacegroup: SpaceGroup | MagneticSpaceGroup):
    #     """ Set the spacegroup"""
    #     self._spacegroup = spacegroup
    #     self._build_sites()

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

    def exchange_constraints(self,
                             site_1: LatticeSite | str | ArrayLike,
                             site_2: LatticeSite | str | ArrayLike):
        """ Get the constraints """
        if isinstance(site_1, str):
            site_1 = self.site_by_name(site_1)

        if isinstance(site_2, str):
            site_2 = self.site_by_name(site_2)

        if not isinstance(site_1, LatticeSite):
            try:
                site_1 = LatticeSite(i=float(site_1[0]),
                                     j=float(site_1[1]),
                                     k=float(site_1[2]),
                                     name="tmp_site_1")
            except Exception as e:
                raise TypeError("Expected `site_1` to be a LatticeSite, vector or a name") from e

        if not isinstance(site_2, LatticeSite):
            try:
                site_2 = LatticeSite(i=float(site_2[0]),
                                     j=float(site_2[1]),
                                     k=float(site_2[2]),
                                     name="tmp_site_1")
            except Exception as e:
                raise TypeError("Expected `site_2` to be a LatticeSite, vector or a name") from e

        return self.spacegroup.exchange_constraints(site_1, site_2)



