"""Space groups and magnetic space groups."""
from collections import defaultdict

from abc import ABC, abstractmethod
from typing import Iterator

import numpy as np
import spglib
from difflib import get_close_matches

from numpy._typing import ArrayLike

from pyspinw.cell_offsets import CellOffsetCoercible, CellOffset
from pyspinw.serialisation import SPWSerialisable, SPWSerialisationContext, SPWDeserialisationContext, expects_keys
from pyspinw.site import Site, ImpliedLatticeSite
from pyspinw.symmetry.canonise import canonise_string
from pyspinw.symmetry.spacegroup_lookup import canonical_aliases, canonical_to_formatted, preferred_names

from pyspinw.symmetry.operations import MagneticOperation, SpaceOperation
from pyspinw.symmetry.settings import Setting
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.symmetry.supercell import Supercell
from pyspinw.symmetry.system import LatticeSystem, lattice_system_letter_lookup, \
    lattice_system_name_lookup
from pyspinw.symmetry.symmetry_checking import ExchangeMatrixConstraints, AnisotropyMatrixConstraints

from pyspinw.tolerances import tolerances

class SymmetryGroup(ABC):
    """Base class for symmetry group and magnetic symmetry group."""

    @property
    @abstractmethod
    def space_operations(self) -> list[SpaceOperation]:
        """The SpaceOperations for this group."""


    def implied_sites_for(self, site: Site, unit_cell: UnitCell) -> list[ImpliedLatticeSite]:
        """Find "duplicate" sites of a given site.

        Parameters
        ----------
        site : Site
            Parent site to generate symmetry-equivalent sites from.
        unit_cell : UnitCell
            Unit cell defining the coordinate system for spin transformations.

        Returns
        -------
        list of ImpliedLatticeSite
            List of symmetry-implied lattice sites distinct from `site`.
        """
        coordinates = site.ijk.reshape(1, 3) % 1
        spins = unit_cell.spin_cartesian_to_lattice_units(site.spin_data)

        new_coordinates = []
        for operation in self.operations:

            candidate, candidate_spins = operation.transform_positions_and_spins(coordinates, spins)

            # If its not the input, continue
            if np.all(np.abs(candidate - coordinates) < tolerances.SAME_SITE_ABS_TOL):
                continue

            # Is it one we've already found
            new = True
            for ijk, spin in new_coordinates:
                if np.all(np.abs(candidate - ijk) < tolerances.SAME_SITE_ABS_TOL):
                    new = False
                    break

            if new:
                new_coordinates.append((candidate, candidate_spins))


        new_sites = []
        for i, (coordinates, spin_data) in enumerate(new_coordinates):
            new_site = ImpliedLatticeSite(
                i=coordinates[0,0],
                j=coordinates[0,1],
                k=coordinates[0,2],
                supercell_spins=unit_cell.spin_lattice_units_to_cartesian(spin_data),
                parent_site=site,
                name=site.name + f" [{i+1}]",
                metadata=site.metadata)

            new_sites.append(new_site)

        return new_sites

    def operations_between_sites(self,
                                 site_1: "Site",
                                 site_2: "Site",
                                 offset: tuple[int, int, int] = (0,0,0),
                                 tolerance=1e-10) -> list[SpaceOperation]:
        """Get a list of symmetry operations that can transform `site_1` into `site_2`.

        Parameters
        ----------
        site_1 : Site
            Starting site.
        site_2 : Site
            Target site.
        offset : tuple of (int, int, int) or CellOffset, optional
            Unit cell offset added to `site_2`. Default is ``(0, 0, 0)``.
        tolerance : float, optional
            Absolute tolerance for comparing site positions. Default is ``1e-10``.

        Returns
        -------
        list of SpaceOperation
            Symmetry operations mapping `site_1` to `site_2` plus `offset`.
        """
        offset = CellOffset.coerce(offset)

        return [operation for operation in self.space_operations
                if np.allclose(operation.transform_positions([site_1.ijk]),
                               [site_2.ijk + offset.vector],
                               atol=tolerance)]

    def operations_on_single_site_pairs(self,
                                        site_1: Site,
                                        site_2: Site,
                                        offset: CellOffsetCoercible = (0, 0, 0)) \
            -> tuple[set[SpaceOperation], set[SpaceOperation]]:
        r"""Get operations on pairs of sites.

        Parameters
        ----------
        site_1 : Site
            First site in the pair.
        site_2 : Site
            Second site in the pair.
        offset : CellOffsetCoercible, optional
            Unit cell offset for the pair. Default is ``(0, 0, 0)``.

        Returns
        -------
        identity_operations : set of SpaceOperation
            Operations that preserve both sites in place ($J = M J M^T$).
        inversion_operations : set of SpaceOperation
            Operations that exchange `site_1` and `site_2` ($J = M J^T M^T$).
        """
        offset = CellOffset.coerce(offset)

        #
        # There are two cases here:
        #  1) where the "bond" is preserved by keeping the sites the same
        #  2) where it is preserved but with exchanging the sites
        #
        # In case 1, the constraint is on J, in the form J = M J M^T
        # In case 2, the constraint is on J^T, as J->J^T constitutes an inversion, J = M J^T M^T

        # Case 1:
        identity_operations = set(self.operations_between_sites(site_1, site_1, offset)).intersection(
            set(self.operations_between_sites(site_2, site_2, offset)))

        # Case 2:
        inversion_operations = set(self.operations_between_sites(site_1, site_2, offset)).intersection(
            set(self.operations_between_sites(site_2, site_1, offset)))

        return identity_operations, inversion_operations

    def operations_between_pairs(self,
                                 pair_1: tuple[Site, Site],
                                 pair_2: tuple[Site, Site],
                                 offset: CellOffsetCoercible = (0, 0, 0),
                                 tolerance: float=1e-10) -> set[SpaceOperation]:
        """Operations in this group that transform one ordered pair into another.

        Parameters
        ----------
        pair_1 : tuple of (Site, Site)
            Starting ordered pair of sites.
        pair_2 : tuple of (Site, Site)
            Target ordered pair of sites.
        offset : CellOffsetCoercible, optional
            Unit cell offset applied to `pair_2`. Default is ``(0, 0, 0)``.
        tolerance : float, optional
            Absolute tolerance for position matching. Default is ``1e-10``.

        Returns
        -------
        set of SpaceOperation
            Set of symmetry operations mapping `pair_1` onto `pair_2`.
        """
        offset = CellOffset.coerce(offset)

        left_operations = self.operations_between_sites(pair_1[0], pair_2[0], offset, tolerance=tolerance)
        right_operations = self.operations_between_sites(pair_1[1], pair_2[1], offset, tolerance=tolerance)

        return set(left_operations).intersection(right_operations)




    def exchange_constraints(self,
                             site_1: Site | ArrayLike,
                             site_2: Site | ArrayLike,
                             unit_cell: UnitCell,
                             do_print: bool = True) -> ExchangeMatrixConstraints:
        """Get the details of the allowed exchange matrices.

        Parameters
        ----------
        site_1 : Site or array_like
            First site or 3-vector of fractional coordinates.
        site_2 : Site or array_like
            Second site or 3-vector of fractional coordinates.
        unit_cell : UnitCell
            Unit cell used to convert operations to Cartesian coordinates.
        do_print : bool, optional
            If ``True``, print a summary of the constraints. Default is ``True``.

        Returns
        -------
        ExchangeMatrixConstraints
            Object describing the symmetry constraints on the exchange matrix.

        Raises
        ------
        TypeError
            If `site_1` or `site_2` cannot be converted to a ``Site``.
        """
        if not isinstance(site_1, Site):
            try:
                site_1 = Site(i=float(site_1[0]),
                              j=float(site_1[1]),
                              k=float(site_1[2]),
                              name="tmp_site_1")
            except Exception as e:
                raise TypeError("Expected `site_1` to be a LatticeSite or vector") from e

        if not isinstance(site_2, Site):
            try:
                site_2 = Site(i=float(site_2[0]),
                              j=float(site_2[1]),
                              k=float(site_2[2]),
                              name="tmp_site_1")
            except Exception as e:
                raise TypeError("Expected `site_2` to be a LatticeSite or vector") from e

        identity_operations, inversion_operations = self.operations_on_single_site_pairs(site_1, site_2)

        identity_operations = [op for op in identity_operations if op.symmorphic]
        inversion_operations = [op for op in inversion_operations if op.symmorphic]

        # # For debugging
        # print("\n\n\n")
        # print("Identity operations:", ", ".join(["{%s}"%op for op in identity_operations]))
        # print("Inversion operations:", ", ".join(["{%s}"%op for op in inversion_operations]))

        inversion_transforms = [op.point_operation_in_cartesian(unit_cell) for op in inversion_operations]
        identity_transforms = [op.point_operation_in_cartesian(unit_cell) for op in identity_operations]

        check = ExchangeMatrixConstraints(identity_transforms, inversion_transforms)

        if do_print:
            check.print_summary()

        return check

    def anisotropy_constraints(self,
                               site: Site | ArrayLike,
                               unit_cell: UnitCell,
                               do_print: bool = True) -> AnisotropyMatrixConstraints:
        """Get the details of the allowed anisotropy matrices.

        Parameters
        ----------
        site : Site or array_like
            Site or 3-vector of fractional coordinates.
        unit_cell : UnitCell
            Unit cell used to convert operations to Cartesian coordinates.
        do_print : bool, optional
            If ``True``, print a summary of the constraints. Default is ``True``.

        Returns
        -------
        AnisotropyMatrixConstraints
            Object describing the symmetry constraints on the anisotropy matrix.

        Raises
        ------
        TypeError
            If `site` cannot be converted to a ``Site``.
        """
        if not isinstance(site, Site):
            try:
                site = Site(i=float(site[0]),
                            j=float(site[1]),
                            k=float(site[2]),
                            name="tmp_site")
            except Exception as e:
                raise TypeError("Expected `site_1` to be a LatticeSite or vector") from e

        operations = self.operations_between_sites(site, site)
        operations = [op for op in operations if op.symmorphic]

        # # For debugging
        # print("\n\n\n")
        # print("Operations:", ", ".join(["{%s}"%op for op in operations]))

        transforms = [op.point_operation_in_cartesian(unit_cell) for op in operations]

        check = AnisotropyMatrixConstraints(transforms)

        if do_print:
            check.print_summary()

        return check

    @property
    def lattice_system(self):
        """Get the lattice system object."""
        raise NotImplementedError("`lattice_system` is not implemented in base class")

    def _description(self):
        """Description string for this spacegroup.

        Raises
        ------
        NotImplementedError
            Always raised in the abstract base class.
        """
        raise NotImplementedError("`description()` is not implemented in base class")

    def for_supercell(self, supercell: Supercell):
        """Get the symmetry group of a supercell, as implied by the symmetry of the unit cell.

        Parameters
        ----------
        supercell : Supercell
            Supercell for which to determine the implied symmetry group.

        Returns
        -------
        SpaceGroup
            Implied space group of the supercell.
        """
        return database.spacegroup_by_name("p1") # TODO: Implement properly


class SpaceGroup(SymmetryGroup, SPWSerialisable):
    """Representation of a space group.

    Parameters
    ----------
    hall_number : int
        Hall number (1 to 530) identifying the specific space group setting.
    number : int
        International Tables space group number (1 to 230).
    international_symbol : str
        Full international Hermann-Mauguin symbol.
    short_symbol : str
        Short international Hermann-Mauguin symbol.
    preferred_symbol : str
        Preferred symbol used for serialisation and display.
    operations : list of SpaceOperation
        Symmetry operations belonging to this space group.
    lattice_system : LatticeSystem
        Crystal lattice system associated with this space group.
    choice : str or None
        Setting choice string (e.g. origin or axis choice), or ``None``.
    """

    serialisation_name = "space_group"

    def __init__(self,
                 hall_number,
                 number,
                 international_symbol,
                 short_symbol,
                 preferred_symbol,
                 operations,
                 lattice_system: LatticeSystem,
                 choice: str | None):

        self.hall_number = hall_number
        self.number = number
        self.symbol = international_symbol
        self.short_symbol = short_symbol
        self.preferred_symbol = preferred_symbol
        self.operations = operations
        self._lattice_system = lattice_system
        self.choice = choice
        self.setting = Setting.from_optional_string(choice)

        # This is slightly unusual, make a reference to the lattice system create_unit_cell method
        self.create_unit_cell = lattice_system.create_unit_cell

    @property
    def lattice_system(self):
        """Get the lattice system associated with this spacegroup."""
        return self._lattice_system

    @property
    def space_operations(self):
        """The SpaceOperations for this group."""
        return self.operations


    def _serialisation_string(self):
        """Name to use to refer to this group in serialisation.

        Returns
        -------
        str
            Preferred symbol string for serialising this space group.
        """
        return self.preferred_symbol

    def _serialise(self, context: SPWSerialisationContext):
        return {"name": self._serialisation_string()}

    @staticmethod
    @expects_keys("name")
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        return database.spacegroup_by_name(json["name"])

    @property
    def name(self) -> str:
        """Just the space group name with setting choice if relevant."""
        if self.choice is None:
            return self.symbol
        else:
            return f"{self.symbol} [{self.choice}]"

    def _description(self):
        """Description string for this spacegroup."""
        return self.name

    def __repr__(self):
        if self.choice is None:
            return f"SpaceGroup({self.number}, {self.symbol}, hall={self.hall_number})"
        else:
            return f"SpaceGroup({self.number}, {self.symbol} [{self.choice}], hall={self.hall_number})"

    def __contains__(self, item: SpaceOperation):
        h = hash(item)
        for operation in self.operations:
            if hash(operation) == h:
                return True
        return False

    def __iter__(self) -> Iterator[SpaceOperation]:
        return iter(self.operations)

#
#
# Database stuff
#
#

class ExactMatch:
    """Result for exact matches in database searches.

    Parameters
    ----------
    spacegroup : SpaceGroup
        The exactly matched space group.
    """

    def __init__(self, spacegroup: SpaceGroup):
        self.spacegroup = spacegroup

class PartialMatch:
    """Result for database searches that are not exact.

    Parameters
    ----------
    spacegroups : list of SpaceGroup
        List of partially matching space groups.
    """

    def __init__(self, spacegroups: list[SpaceGroup]):
        self.spacegroups = spacegroups


class NoSuchGroup(Exception):
    """Raised when a group is not found."""


class SpacegroupDatabase:
    """Holds all the information about spacegroups and magnetic spacegroups."""

    # Letter representing the crystal system monoclinic, cubic etc
    _spacegroup_number_to_lattice_system = \
        ["a" for _ in range(1, 3)] + \
        ["m" for _ in range(3, 16)] + \
        ["o" for _ in range(16, 75)] + \
        ["t" for _ in range(75, 143)] + \
        ["h" for _ in range(143, 195)] + \
        ["c" for _ in range(195, 231)]

    # Converts first letter of spacegroup to the relevant Bravais notation
    _spacegroup_symbol_to_bravais_symbol = {
        "A": "S",
        "B": "S",
        "C": "S",
        "F": "F",
        "I": "I",
        "P": "P",
        "R": "R",
    }

    def __init__(self):

        # Create spacegroups
        self._lattice_symbol_to_spacegroups = defaultdict(list[SpaceGroup])
        self.spacegroups: list[SpaceGroup] = []
        self._operation_lookup = []

        for hall_number in range(1, 531):

            # Get the relevant data

            sg_data = spglib.get_spacegroup_type(hall_number)
            op_data = spglib.get_symmetry_from_database(hall_number)

            # Classify
            name = sg_data.international_full
            short_name = sg_data.international_short
            lattice_type_from_name = name[0]

            first_letter = SpacegroupDatabase._spacegroup_number_to_lattice_system[sg_data.number-1]
            second_letter = SpacegroupDatabase._spacegroup_symbol_to_bravais_symbol[lattice_type_from_name]
            bravais_lattice = first_letter + second_letter

            # Load the operations
            translations = op_data["translations"]
            rotations = op_data["rotations"]

            operations = []
            for translation, rotation, in zip(translations, rotations):
                op = SpaceOperation.from_numpy(rotation.T, translation)
                operations.append(op)

            self._operation_lookup.append(
                {canonise_string(operation.text_form) for operation in operations})

            choice = None if sg_data.choice == "" else sg_data.choice

            if choice == "R":
                lattice_system = lattice_system_name_lookup["Rhombohedral"]
                bravais_lattice = "hR"
            else:
                if first_letter == "m":
                    if choice.startswith("a") or choice.startswith("-a"):
                        lattice_system = lattice_system_name_lookup["MonoclinicA"]
                    elif choice.startswith("c") or choice.startswith("-c"):
                        lattice_system = lattice_system_name_lookup["MonoclinicC"]
                    else:
                        lattice_system = lattice_system_name_lookup["MonoclinicB"]

                else:
                    lattice_system = lattice_system_letter_lookup[first_letter]

            group = SpaceGroup(
                hall_number=hall_number,
                number=sg_data.number,
                international_symbol=name,
                short_symbol=short_name,
                preferred_symbol=preferred_names[hall_number],
                operations=operations,
                choice=choice,
                lattice_system=lattice_system)

            self.spacegroups.append(group)
            self._lattice_symbol_to_spacegroups[bravais_lattice].append(group)

        # Lookup for spacegroups by long/short international symbol
        self._spacegroup_symbol_lookup: dict[str, SpaceGroup] = {group.symbol: group for group in self.spacegroups}
        self._spacegroup_symbol_lookup.update({group.short_symbol: group for group in self.spacegroups})

        # Lookup for spacegroups based on their number
        self._spacegroup_number_lookup = defaultdict(list[SpaceGroup])
        for group in self.spacegroups:
            self._spacegroup_number_lookup[group.number].append(group)



    def all_spacegroup_settings(self, number: int):
        """Get a list of all spacegroups that correspond to a given setting.

        Parameters
        ----------
        number : int
            International Tables space group number (1 to 230).

        Returns
        -------
        list of SpaceGroup
            All space group settings corresponding to `number`.
        """
        return self._spacegroup_number_lookup[number]

    def spacegroups_for_lattice_symbol(self, lattice_symbol: str):
        """Get the spacegroups that correspond to a given lattice symbol (e.g. hR).

        Parameters
        ----------
        lattice_symbol : str
            Two-letter Pearson/Bravais lattice symbol (e.g. ``'hR'``).

        Returns
        -------
        list of SpaceGroup
            Space groups associated with `lattice_symbol`.
        """
        return self._lattice_symbol_to_spacegroups[lattice_symbol]

    def spacegroups_with_operations(self, operation_string: str) -> ExactMatch | PartialMatch:
        """Get a spacegroup by a list of operations.

        Example: ``"x, y, z; -x, -y, -z"`` should give P1.

        Parameters
        ----------
        operation_string : str
            A list of operations in the ``x, y, z`` form, semicolon separated.

        Returns
        -------
        ExactMatch or PartialMatch
            ``ExactMatch`` if a group has the exact set of operations, or
            ``PartialMatch`` containing groups that are supersets of the operations.
        """
        operations = {canonise_string(op) for op in operation_string.split(";")}

        out = []
        for i, test_operations in enumerate(self._operation_lookup):
            if operations <= test_operations:
                if operations == test_operations:
                    return ExactMatch(self.spacegroups[i])

                out.append(self.spacegroups[i])

        return PartialMatch(out)


    def spacegroup_by_name(self, name: str) -> SpaceGroup:
        """Get a spacegroup by name.

        Parameters
        ----------
        name : str
            Name or alias of the space group.

        Returns
        -------
        SpaceGroup
            Matching ``SpaceGroup`` instance.

        Raises
        ------
        NoSuchGroup
            If no space group matches `name`.
        """
        canonised_input = canonise_string(name)

        if canonised_input in canonical_aliases:
            index = canonical_aliases[canonised_input]
            return self.spacegroups[index-1]

        # Failure branch...
        # We want to get names that are similar in their input text form, but report the standard form
        similar = get_close_matches(name, canonical_aliases.keys(), n=3, cutoff=0.4)

        if similar:

            #

            suggestion_string = ", ".join([f"'{canonical_to_formatted[s]}'" for s in similar])
            message_string = (f"Unknown space group '{name}', "
                              f"perhaps you meant {suggestion_string} or something similar.")

        else:
            message_string = f"Unknown space group '{name}', doesn't seem to be even close to a magnetic space group."

        raise NoSuchGroup(message_string)

# Create the database
database = SpacegroupDatabase()


if __name__ == "__main__":
    # print(database.spacegroup_by_name("b2/m"))
    print(database.spacegroup_by_name("c2/m"))

    print(database.spacegroup_by_name("pnnn:1"))

    for i, spacegroup in enumerate(database.spacegroups):
        print(i, spacegroup)
