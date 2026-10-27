"""Magnetic structures."""
import logging

import numpy as np
from ase.data import chemical_symbols
from numpy._typing import ArrayLike

from pyspinw.cell_offsets import CellOffset
from pyspinw.exchangegroup import DirectionalityFilter
from pyspinw.lattice_distances import full_search_space
from pyspinw.cell_offsets import cell_offset_generator
from pyspinw.serialisation import SPWSerialisable, SPWSerialisationContext, SPWDeserialisationContext, expects_keys
from pyspinw.site import Site
from pyspinw.symmetry.group import SpaceGroup, SymmetryGroup, database
from pyspinw.symmetry.operations import SpaceOperation
from pyspinw.symmetry.supercell import Supercell, TiledSupercell, DirectSupercell
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.tolerances import tolerances
from pyspinw.util import connected_components, arraylike_equality, IncrementalPointHistogram

logger = logging.getLogger("Structure")

class Structure(SPWSerialisable):
    """Representation of the magnetic structure.

    Parameters
    ----------
    sites : list of Site
        Sites used to define the structure.
    unit_cell : UnitCell
        Unit cell of the crystal structure.
    spacegroup : SymmetryGroup or None, optional
        Space group or symmetry group of the structure. If ``None``, P1
        symmetry is used. Default is ``None``.
    supercell : Supercell or None, optional
        Supercell definition. If ``None``, a ``TiledSupercell`` is used.
        Default is ``None``.
    skip_checks : bool, optional
        Whether to skip validation checks on supercell spin dimensions and unit
        cell symmetry constraints. Default is ``False``.
    show_unit_cell_warning : bool, optional
        Whether to log a warning if the unit cell violates negative lattice
        constraints. Default is ``True``.
    assume_symmetric_sites_are_provided : bool, optional
        Whether ``symmetric_sites`` already contains all symmetry-implied
        sites. Default is ``False``.
    cooerce_spin_data_to_match_supercell : bool, optional
        Whether single-component site spin data should be repeated to match
        the number of supercell components. Default is ``True``.
    symmetric_sites : list of Site or None, optional
        Pre-computed list of symmetry-expanded sites used when
        ``assume_symmetric_sites_are_provided`` is ``True``. Default is
        ``None``.
    """

    serialisation_name = "structure"

    def __init__(self,
                 sites: list[Site],
                 unit_cell: UnitCell,
                 spacegroup: SymmetryGroup | None = None,
                 supercell: Supercell | None = None,
                 skip_checks: bool = False,
                 show_unit_cell_warning: bool=True,
                 assume_symmetric_sites_are_provided: bool=False,
                 cooerce_spin_data_to_match_supercell: bool=True,
                 symmetric_sites: list[Site] | None = None):


        if not isinstance(unit_cell, UnitCell):
            raise TypeError("Expected `unit_cell` to be of type `UnitCell`")

        if spacegroup is not None:
            if not isinstance(spacegroup, SymmetryGroup):
                raise TypeError("Expected `spacegroup` to be of type `SymmetryGroup` (or SpaceGroup)")

        if supercell is not None:
            if not isinstance(supercell, Supercell):
                raise TypeError("Expected `supercell` to be a subclass of `Supercell`")

        self._spacegroup = database.spacegroups[0] if spacegroup is None else spacegroup

        self._input_sites = sites
        self._input_uid_to_site = {site.unique_id: site for site in sites}
        self._unit_cell = unit_cell

        self._supercell = TiledSupercell() if supercell is None else supercell

        if assume_symmetric_sites_are_provided and symmetric_sites is not None:
            self._sites = symmetric_sites
        else:
            self._sites: list[Site] = self._extended_sites()

        self._site_lookup = {site.unique_id: site for site in self._sites} # Cache a map

        if not skip_checks:

            # Check that supercell components match site dimensions
            n_components = self._supercell.n_components()
            maybe_warn = False
            bad_sites = []
            for site in self.sites:
                if site.n_components() != n_components:

                    if cooerce_spin_data_to_match_supercell and site.n_components() == 1:
                        site.spin_data = np.repeat(site.spin_data, n_components, axis=0)
                        maybe_warn = True
                    else:
                        bad_sites.append(site)

            if maybe_warn and not isinstance(self._supercell, DirectSupercell):
                logger.warning("`spin_data` for some sites has been adjusted to match the supercell")

            if bad_sites:
                raise ValueError("Expected the shape of site spin data to match what the supercell requires "
                                 f"({supercell.n_components()}-by-3), "
                                 "bad sites are: " + ", ".join([
                        site.name if site.name is not None and site.name != "" else f"UID-{site.unique_id}"
                          for site in bad_sites]))

            # check that the unit cell is consistent with the spacegroup

            if self._spacegroup.lattice_system.constrain(unit_cell) != unit_cell:
                raise ValueError(f"{unit_cell} is not compatible with spacegroup {self._spacegroup}, "
                                 f"requires lattice type {self._spacegroup.lattice_system.name}")

            # Check that the unit cell fulfils all the constraints of the unit cell
            if show_unit_cell_warning and spacegroup is not None:
                violations = self._spacegroup.lattice_system.violated_negative_constraints(unit_cell)
                if violations:
                    logger.warning(f"{unit_cell} violates the conditions {violations} for "
                                   f"{self._spacegroup.lattice_system} used by {self._spacegroup}")



    def full_structure_site_list(self):
        """All the sites in the structure.

        Returns
        -------
        list of Site
            All sites across all cells in the supercell, including
            symmetry-implied sites.
        """
        all_sites = []

        # unit cell stuff
        cell_sites = self._extended_sites()

        # supercell
        for cell in self.supercell.cells():
            all_sites += [Site(
                i=site.i + cell.i,
                j=site.j + cell.j,
                k=site.k + cell.k,
                supercell_spins=site.spin_data,
                name=site.name
            ) for site in cell_sites]

        return all_sites

    def matplotlib_site_data(self):
        """Data for making matplotlib scatter plots of sites.

        Returns
        -------
        x : numpy.ndarray
            Cartesian x-coordinates of all sites in the supercell.
        y : numpy.ndarray
            Cartesian y-coordinates of all sites in the supercell.
        z : numpy.ndarray
            Cartesian z-coordinates of all sites in the supercell.
        """
        array_data = np.array([self.unit_cell.lattice_units_to_cartesian(site._ijk)
                         for site in self.full_structure_site_list()])

        return array_data[:,0], array_data[:, 1], array_data[:, 2]

    def _extended_sites(self)  -> list[Site]:
        """All the sites, including those implied by symmetry.

        Returns
        -------
        list of Site
            Unique sites within the unit cell after applying space group
            operations and resolving site collisions.

        Raises
        ------
        ValueError
            If multiple input sites occupy the same spatial location.
        """
        site_list = self._input_sites.copy()
        for site in self._input_sites:
            site_list += self._spacegroup.implied_sites_for(site, self.unit_cell)

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
                unique_sites.append(Site(
                    site_1.i,
                    site_1.j,
                    site_1.k,
                    supercell_spins=site_1.spin_data,
                    name = site_1.name,
                    metadata = site_1.metadata.copy()
                ))

            else:
                unique_sites.append(Site(
                    site_1.i,
                    site_1.j,
                    site_1.k,
                    supercell_spins=np.zeros_like(site_1.spin_data),
                    name=site_1.name,
                    metadata=site_1.metadata.copy()
                ))

        return unique_sites

    def _expansion_site_mapping(self, supercell_size: tuple[int, int, int] | None = None):
        """Expand supercell into a single, bigger cell.

        Parameters
        ----------
        supercell_size : tuple of (int, int, int) or None, optional
            Explicit supercell scaling ``(na, nb, nc)``. If ``None``, the
            current supercell's ``cell_size()`` is used. Default is ``None``.

        Returns
        -------
        big_cell : UnitCell
            Expanded unit cell scaled by ``supercell_size``.
        mapping : dict
            Mapping from ``(site_unique_id, cell_offset_tuple)`` to the
            corresponding ``Site`` in the expanded unit cell.
        """
        # Calculate new cell
        if supercell_size is None:
            scale = self.supercell.cell_size()
        else:
            scale = supercell_size

        big_cell = self.unit_cell.updated(
            a=self.unit_cell.a * scale[0],
            b=self.unit_cell.b * scale[1],
            c=self.unit_cell.c * scale[2])

        scaling = np.array(scale, dtype=float)

        # Create a mapping between sites and offsets to the new sites
        mapping: dict[tuple[int, tuple[int, int, int]], Site] = {}
        for index, offset in enumerate(cell_offset_generator(*scale)):
            for site in self.sites:

                position = offset.vector + np.array(site.ijk, dtype=float)
                position /= scaling


                spin = self.supercell.spin(site, cell_offset=offset)

                new_site = Site(
                    i=position[0],
                    j=position[1],
                    k=position[2],
                    supercell_spins=spin,
                    g=site.g,
                    name=f"{site.name}[{index}]",
                    metadata=site.metadata.copy())

                mapping[(site._unique_id, offset.as_tuple)] = new_site

        return big_cell, mapping

    def site_pairs(self):
        """Generator for pairs of sites.

        Yields
        ------
        site_1 : Site
            First site in the pair.
        site_2 : Site
            Second site in the pair (distinct from ``site_1``).
        """
        for index_1, site_1 in enumerate(self.sites):
            for index_2, site_2 in enumerate(self.sites):
                if index_1 == index_2:
                    continue
                yield site_1, site_2

    def expand(self):
        """Expand supercell into a single, bigger cell.

        Returns
        -------
        Structure
            New structure whose unit cell encompasses the entire supercell.
        """
        cell, mapping = self._expansion_site_mapping()


        return Structure(
            sites=[site for site in mapping.values()],
            unit_cell=cell,
            spacegroup=self.spacegroup.for_supercell(self.supercell),
            supercell=TiledSupercell(scaling=(1, 1, 1))
        )

    def _without_nonmagnetic_with_removed_uids(self) -> tuple["Structure", set[int]]:
        """Get a copy of the structure without non-magnetic sites, along with a list of uids for the removed sites.

        Returns
        -------
        new_structure : Structure
            Structure with non-magnetic sites removed.
        excluded_uids : set of int
            Unique IDs of the removed non-magnetic sites.
        """
        new_sites: list[Site] = []

        # Get rid of the nonmagnetic atoms from input sites
        new_input_sites: list[Site] = []
        for site in self._input_sites:
            if not np.allclose(site.spin_data, 0.0):
                new_input_sites.append(site)

        # Get rid of them from all the sites
        excluded_uids: set[int] = set()
        for site in self._sites:
            if np.allclose(site.spin_data, 0.0):
                excluded_uids.add(site.unique_id)
            else:
                new_sites.append(site)

        new_structure = Structure(
            new_input_sites,
            self._unit_cell,
            self._spacegroup,
            self._supercell,
            show_unit_cell_warning = False,
            assume_symmetric_sites_are_provided=True,
            symmetric_sites = new_sites
            )

        return new_structure, excluded_uids

    def without_nonmagnetic(self) -> "Structure":
        """Get a copy of this structure with non-magnetic sites removed.

        Returns
        -------
        Structure
            Copy of this structure containing only magnetic sites.
        """
        structure, _ = self._without_nonmagnetic_with_removed_uids()

        return structure

    def _build_sites(self):
        """Updates the site variable on parameter changes."""
        self._sites = self._extended_sites()

    @property
    def sites(self) -> list[Site]:
        """Get the sites used to define the structure and implied by symmetry."""
        return self._sites.copy()
        # return self._input_sites.copy()

    def sites_by_name(self, name) -> list[Site]:
        """Find sites by their name, can make multiple matches.

        Parameters
        ----------
        name : str
            Case-insensitive prefix used to match site names.

        Returns
        -------
        list of Site
            Sites whose names start with ``name`` (case-insensitive).
        """
        return [site for site in self._sites if site.name.lower().startswith(name.lower())]

    def site_by_name(self, name: str):
        """Get a single site by its name.

        Parameters
        ----------
        name : str
            Site name or prefix to look up.

        Returns
        -------
        Site
            Matching site.

        Raises
        ------
        ValueError
            If no site matches ``name`` or if multiple sites match without a
            unique exact match.
        """
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

    def sites_by_element(self, element: str | None) -> list[Site]:
        """Get list of sites with the specified element in the metadata.

        Parameters
        ----------
        element : str or None
            Chemical element symbol to match (or ``None``).

        Returns
        -------
        list of Site
            Sites whose metadata element equals ``element``.

        Raises
        ------
        ValueError
            If ``element`` is not a known chemical symbol.
        """
        if element is not None and element not in chemical_symbols[1:]:
            raise ValueError(f"{element} is not a known element")
        return [site for site in self._sites if site.metadata.element == element]

    def without_elements(self, *elements):
        """Get the structure without sites with given elements.

        Parameters
        ----------
        *elements : str or None
            Chemical element symbols to exclude from the structure.

        Returns
        -------
        Structure
            New structure with sites matching ``elements`` removed.

        Raises
        ------
        ValueError
            If any specified element is not a known chemical symbol.
        """
        for element in elements:
            if element is not None and element not in chemical_symbols[1:]:
                raise ValueError(f"{element} is not a known element")

        new_sites = [site for site in self._sites if site.metadata.element not in elements]

        return Structure(new_sites, unit_cell=self.unit_cell, spacegroup=self.spacegroup, supercell=self.supercell)


    @property
    def text_summary(self) -> str:
        """Textual details of this structure."""
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
        """Print out details of this structure."""
        print(self.text_summary)

    def all_neighbours(self,
                       site: Site,
                       n=1,
                       element: str | None = None,
                       direction_filter: DirectionalityFilter | None = None) \
            -> list[list[tuple[Site, CellOffset]]]:
        """Get list-of-lists containing (i<=n)-th nearest neighbours along with their cell offsets.

        Parameters
        ----------
        site : Site
            Reference site from which neighbour distances are measured.
        n : int, optional
            Maximum neighbour order ``n`` for ``n``-th nearest neighbours.
            Default is ``1``.
        element : str or None, optional
            Restrict search to sites of this chemical element. Default is
            ``None``.
        direction_filter : DirectionalityFilter or None, optional
            Filter the results using this directionality filter. Default is
            ``None``.

        Returns
        -------
        list of list of tuple of (Site, CellOffset)
            List of ``(site, offset)`` pairs for each neighbour order from 0
            up to and including ``n``.

        Raises
        ------
        ValueError
            If ``element`` is not a valid chemical element symbol.
        """
        if element is not None and element not in chemical_symbols[1:]:
            raise ValueError(f"{element} is not an element")

        # Get a list of sites we want to check
        sites_to_check: list[Site] = []
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
                   site: Site,
                   n=1,
                   element: str | None = None,
                   direction_filter: DirectionalityFilter | None = None):
        """Get a list of n-th nearest neighbours, along with cell offsets.

        Parameters
        ----------
        site : Site
            Reference site from which neighbour distances are measured.
        n : int, optional
            Neighbour order ``n`` for ``n``-th nearest neighbours. Default is
            ``1``.
        element : str or None, optional
            Restrict search to sites of this chemical element. Default is
            ``None``.
        direction_filter : DirectionalityFilter or None, optional
            Filter the results using this directionality filter. Default is
            ``None``.

        Returns
        -------
        list of tuple of (Site, CellOffset)
            List of ``(site, offset)`` pairs for the ``n``-th nearest
            neighbours.
        """
        all_neighbors = self.all_neighbours(
                site = site,
                n=n,
                element=element,
                direction_filter=direction_filter)

        return all_neighbors[-1]

    def one_neighbour(self,
                      site: Site,
                      n=1,
                      element: str | None = None,
                      direction_filter: DirectionalityFilter | None = None):
        """Get one, arbitrary, n-th nearest neighbour, along with cell offset.

        Parameters
        ----------
        site : Site
            Reference site from which neighbour distances are measured.
        n : int, optional
            Neighbour order ``n`` for ``n``-th nearest neighbours. Default is
            ``1``.
        element : str or None, optional
            Restrict search to sites of this chemical element. Default is
            ``None``.
        direction_filter : DirectionalityFilter or None, optional
            Filter the results using this directionality filter. Default is
            ``None``.

        Returns
        -------
        tuple of (Site, CellOffset)
            Single ``(site, offset)`` pair for an ``n``-th nearest neighbour.

        Raises
        ------
        ValueError
            If no neighbours pass the filter.
        """
        neighbors = self.neighbours(
            site=site,
            n=n,
            element=element,
            direction_filter=direction_filter)

        if len(neighbors) == 0:
            raise ValueError("No neighbours that pass the filter")

        return neighbors[0]


    def symmetry_related(self, site: Site) -> list[tuple[Site, set[SpaceOperation]]]:
        """Get a list of sites related to the specified site by symmetry, including the original site.

        Parameters
        ----------
        site : Site
            Reference site to match under space group operations.

        Returns
        -------
        list of tuple of (Site, set of SpaceOperation)
            Pairs of ``(other_site, operations)`` where ``operations`` is the
            non-empty set of space group operations mapping ``site`` to
            ``other_site``.
        """
        symmetry_related = []
        for other_site in self._sites:

            ops = self.spacegroup.operations_between_sites(site, other_site)
            if len(ops) > 0:
                symmetry_related.append((other_site, set(ops)))

        return symmetry_related

    @property
    def text_summary(self) -> str:
        """Textual details of this structure."""
        lines = []
        lines.append(f"Unit Cell: {self.unit_cell.text_summary}")
        lines.append(f"Spacegroup: {self.spacegroup._description()}")
        supercell_text_data = self.supercell.text_data()
        lines.append(f"Supercell: {supercell_text_data[0]}")
        lines += ["  " + s for s in supercell_text_data[1:]]
        lines.append("Sites:")
        for site in self._sites:

            is_not_input_chr = "" if site.unique_id in self._input_uid_to_site else "* "

            lines.append(f"  {is_not_input_chr}{site}")

        return "\n".join(lines)

    def print_summary(self):
        """Print out details of this structure."""
        print(self.text_summary)

    @property
    def spacegroup(self) -> SpaceGroup:
        """Get the spacegroup."""
        return self._spacegroup

    @property
    def unit_cell(self) -> UnitCell:
        """Get the unit cell."""
        return self._unit_cell

    @unit_cell.setter
    def unit_cell(self, unit_cell: UnitCell):
        self._unit_cell = unit_cell
        self._build_sites()

    @property
    def supercell(self) -> Supercell:
        """Get the supercell."""
        return self._supercell

    @supercell.setter
    def supercell(self, supercell: Supercell):
        self._supercell = supercell
        self._build_sites()

    def updated(self,
                sites: list[Site] | None = None,
                unit_cell: UnitCell | None = None,
                spacegroup: SymmetryGroup | None = None,
                supercell: Supercell | None = None):
        """Make a copy of this structure, but with selected (not None) fields replaced.

        Parameters
        ----------
        sites : list of Site or None, optional
            Replacement site list. If ``None``, the current sites are reused.
            Default is ``None``.
        unit_cell : UnitCell or None, optional
            Replacement unit cell. If ``None``, the current unit cell is
            reused. Default is ``None``.
        spacegroup : SymmetryGroup or None, optional
            Replacement space group. If ``None``, the current space group is
            reused. Default is ``None``.
        supercell : Supercell or None, optional
            Replacement supercell. If ``None``, the current supercell is
            reused. Default is ``None``.

        Returns
        -------
        Structure
            Updated copy of this structure.
        """
        symmetry_considerations_changed = sites is not None and spacegroup is not None
        sites_to_send = None if symmetry_considerations_changed else self.sites

        return Structure(
            self._input_sites if sites is None else sites,
            self._unit_cell if unit_cell is None else unit_cell,
            self._spacegroup if spacegroup is None else spacegroup,
            self._supercell if supercell is None else supercell,
            show_unit_cell_warning= not (spacegroup is None and unit_cell is None),
            assume_symmetric_sites_are_provided=not symmetry_considerations_changed,
            symmetric_sites = sites_to_send
            )

    def exchange_constraints(self,
                             site_1: Site | str | ArrayLike,
                             site_2: Site | str | ArrayLike,
                             do_print=True):
        """Get the constraints for the exchange matrix between two sites.

        Parameters
        ----------
        site_1 : Site or str or array_like
            First site, specified as a ``Site``, site name, or 3-element
            fractional coordinate vector.
        site_2 : Site or str or array_like
            Second site, specified as a ``Site``, site name, or 3-element
            fractional coordinate vector.
        do_print : bool, optional
            Whether to print the constraints to standard output. Default is
            ``True``.

        Returns
        -------
        ParameterizationResult
            Symmetry constraints on the exchange matrix between ``site_1`` and
            ``site_2``.

        Raises
        ------
        TypeError
            If ``site_1`` or ``site_2`` cannot be resolved to a ``Site``.
        """
        if isinstance(site_1, str):
            site_1 = self.site_by_name(site_1)

        if isinstance(site_2, str):
            site_2 = self.site_by_name(site_2)

        if not isinstance(site_1, Site):
            try:
                site_1 = Site(i=float(site_1[0]),
                              j=float(site_1[1]),
                              k=float(site_1[2]),
                              name="tmp_site_1")
            except Exception as e:
                raise TypeError("Expected `site_1` to be a LatticeSite, vector or a name") from e

        if not isinstance(site_2, Site):
            try:
                site_2 = Site(i=float(site_2[0]),
                              j=float(site_2[1]),
                              k=float(site_2[2]),
                              name="tmp_site_1")
            except Exception as e:
                raise TypeError("Expected `site_2` to be a LatticeSite, vector or a name") from e

        return self.spacegroup.exchange_constraints(site_1, site_2, self.unit_cell, do_print=do_print)


    def anisotropy_constraints(self,
                               site: Site | str | ArrayLike,
                               do_print=True):
        """Get the constraints on anisotropies for a given site.

        Parameters
        ----------
        site : Site or str or array_like
            Target site, specified as a ``Site``, site name, or 3-element
            fractional coordinate vector.
        do_print : bool, optional
            Whether to print the constraints to standard output. Default is
            ``True``.

        Returns
        -------
        ParameterizationResult
            Symmetry constraints on the anisotropy matrix at ``site``.

        Raises
        ------
        TypeError
            If ``site`` cannot be resolved to a ``Site``.
        """
        if isinstance(site, str):
            site = self.site_by_name(site)

        if not isinstance(site, Site):
            try:
                site = Site(i=float(site[0]),
                            j=float(site[1]),
                            k=float(site[2]),
                            name="tmp_site_1")
            except Exception as e:
                raise TypeError("Expected `site` to be a LatticeSite, vector or a name") from e

        return self.spacegroup.anisotropy_constraints(site, self.unit_cell, do_print=do_print)

    def _serialise(self, context: SPWSerialisationContext) -> dict:
        return {
            "input_sites": [site._serialise(context) for site in self._input_sites],
            "unit_cell": self.unit_cell._serialise(context),
            "spacegroup": self.spacegroup._serialise(context),
            "supercell": self.supercell._serialise(context),
            "symmetric_sites": [site._serialise(context) for site in self.sites]
            }

    @staticmethod
    @expects_keys("input_sites, unit_cell, spacegroup, supercell, symmetric_sites")
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        return Structure(sites=[Site._deserialise(s, context) for s in json["input_sites"]],
                         unit_cell=UnitCell._deserialise(json["unit_cell"], context),
                         spacegroup=SpaceGroup._deserialise(json["spacegroup"], context),
                         supercell=Supercell._deserialise(json["supercell"], context),
                         skip_checks= True,
                         show_unit_cell_warning=False,
                         assume_symmetric_sites_are_provided=True,
                         cooerce_spin_data_to_match_supercell=False,
                         symmetric_sites=[Site._deserialise(s, context) for s in json["symmetric_sites"]])

    def __repr__(self):

        parts = [f"{len(self.sites)} sites",
                 self.unit_cell.short_name,
                 self.spacegroup.name]

        if not (isinstance(self.supercell, TiledSupercell) and self.supercell.scaling == (1, 1, 1)):
            parts.append(f"{self.supercell.supercell_name} supercell")


        contents = ", ".join(parts)

        return f"Structure({contents})"
