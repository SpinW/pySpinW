"""Groups of exchanges, a neater way of representing multiple, similar exchanges."""

import numpy as np
from numpy._typing import ArrayLike

from pyspinw.exchange import Exchange
from pyspinw.batch_exchanges import batch_exchanges, default_naming_pattern
from pyspinw.symmetry.symmetry_settings import SymmetrySettings
from pyspinw.site import Site
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.tolerances import tolerances


class AbstractExchangeGroup:
    """Base class for exchange groups."""

    def exchanges(self, sites: list[Site], symmetry: SymmetrySettings):
        """Get the exchanges defined by this groups given the sites and symmetry provided.

        Parameters
        ----------
        sites : list of Site
            List of sites in the structure.
        symmetry : SymmetrySettings
            Symmetry settings of the crystal structure.

        Raises
        ------
        NotImplementedError
            Always raised in the abstract base class.
        """
        raise NotImplementedError("exchanges not available in base class")

class DirectionalityFilter:
    """Base class for filtering exchanges by direction."""

    def accept(self, vector: np.ndarray) -> bool:
        """Return True if in an allowed direction.

        Parameters
        ----------
        vector : (3,) numpy.ndarray
            Bond direction vector in Cartesian coordinates.

        Returns
        -------
        bool
            ``True`` if ``vector`` lies in an allowed direction.

        Raises
        ------
        NotImplementedError
            Always raised in the abstract base class.
        """
        raise NotImplementedError("accept() not implemented in base class")

class NoDirectionalityFilter:
    """Direction filter that always accepts."""

    def accept(self, vector: np.ndarray) -> bool:
        """Return True if in an allowed direction (in this case, always).

        Parameters
        ----------
        vector : (3,) numpy.ndarray
            Bond direction vector in Cartesian coordinates.

        Returns
        -------
        bool
            Always ``True``.
        """
        return True

class InDirectionFilter(DirectionalityFilter):
    """Selects vectors in a given direction.

    Parameters
    ----------
    direction : array_like
        Target 3D direction vector.
    max_dev_angle_deg : float, optional
        Maximum allowed angular deviation in degrees. Default is ``0.01``.

    Attributes
    ----------
    direction : (3,) numpy.ndarray
        Unit vector along the target direction.
    in_direction_dev_num : float
        Cosine of the maximum allowed deviation angle.
    """

    def __init__(self, direction: ArrayLike, max_dev_angle_deg: float = 0.01):
        direction = np.array(direction, dtype=float)
        self.direction = direction / np.sqrt(np.sum(direction**2))
        self.in_direction_dev_num = np.cos((np.pi/180) * max_dev_angle_deg)

    def accept(self, vector):
        """DirectionalityFilter implementation: return True if in a similar direction."""
        sq_mag = np.sum(vector**2)
        if sq_mag == 0:
            return False
        else:
            return np.dot(self.direction, vector) / np.sqrt(sq_mag) > self.in_direction_dev_num


class BiDirectionFilter(DirectionalityFilter):
    """Selects vectors in a given direction.

    Parameters
    ----------
    direction : array_like
        Target 3D axis vector (accepts both parallel and antiparallel vectors).
    max_dev_angle_deg : float, optional
        Maximum allowed angular deviation in degrees. Default is ``0.01``.

    Attributes
    ----------
    direction : (3,) numpy.ndarray
        Unit vector along the target axis.
    in_direction_dev_num : float
        Cosine of the maximum allowed deviation angle.
    """

    def __init__(self, direction: ArrayLike, max_dev_angle_deg: float = 0.01):
        direction = np.array(direction, dtype=float)
        self.direction = direction / np.sqrt(np.sum(direction**2))
        self.in_direction_dev_num = np.cos((np.pi/180) * max_dev_angle_deg)

    def accept(self, vector):
        """DirectionalityFilter implementation: return True if in a similar direction."""
        sq_mag = np.sum(vector**2)
        if sq_mag == 0:
            return False
        else:
            return np.abs(np.dot(self.direction, vector) / np.sqrt(sq_mag)) > self.in_direction_dev_num


class InPlaneFilter(DirectionalityFilter):
    """Selects vectors in a given plane (specified by normal).

    Parameters
    ----------
    direction : array_like
        3D normal vector to the target plane.
    max_dev_angle_deg : float, optional
        Maximum allowed angular deviation from the plane in degrees.
        Default is ``0.01``.

    Attributes
    ----------
    direction : (3,) numpy.ndarray
        Unit normal vector of the plane.
    in_plane_dev_num : float
        Sine of the maximum allowed out-of-plane deviation angle.
    """

    def __init__(self, direction: ArrayLike, max_dev_angle_deg: float = 0.01):
        direction = np.array(direction)
        self.direction = direction / np.sqrt(np.sum(direction**2))
        self.in_plane_dev_num = np.sin((np.pi/180) * max_dev_angle_deg)

    def accept(self, vector):
        """DirectionalityFilter implementation: return True if in the plane normal to provided vector."""
        sq_mag = np.sum(vector ** 2)
        if sq_mag == 0:
            return False
        else:
            return np.abs(np.dot(self.direction, vector) / np.sqrt(sq_mag)) < self.in_plane_dev_num

class ExchangeGroup:
    """Class representing the batch creation of exchanges.

    Parameters
    ----------
    name : str
        Name of the exchange group.
    bond : int
        Bond shell index (if ``> 0``, selects only bonds belonging to this shell).
    min_distance : float
        Minimum bond distance to include when ``bond <= 0``.
    max_distance : float
        Maximum bond distance to include when ``bond <= 0``.
    max_order : int or None
        Maximum shell order to include, or ``None`` for no limit.
    naming_pattern : str or None
        Naming pattern string for generated exchanges, or ``None`` to use
        ``default_naming_pattern``.
    exchange_type : type of Exchange
        `Exchange` subclass to instantiate for each bond.
    exchange_parameters : dict
        Keyword arguments forwarded to ``exchange_type`` upon instantiation.
    direction_filter : DirectionalityFilter or None
        Optional filter to select bonds by direction.

    Attributes
    ----------
    name : str
        Name of the exchange group.
    bond : int
        Bond shell index.
    min_distance : float
        Minimum bond distance.
    max_distance : float
        Maximum bond distance.
    max_order : int or None
        Maximum shell order.
    naming_pattern : str or None
        Naming pattern string.
    exchange_type : type of Exchange
        `Exchange` subclass.
    exchange_parameters : dict
        Keyword arguments for ``exchange_type``.
    direction_filter : DirectionalityFilter or None
        Filter for bond directions.
    """

    def __init__(self,
                 name: str,
                 bond: int,
                 min_distance: float,
                 max_distance: float,
                 max_order: int | None,
                 naming_pattern: str | None,
                 exchange_type: type[Exchange],
                 exchange_parameters: dict,
                 direction_filter: DirectionalityFilter | None):

        self.name = name
        self.bond = bond
        self.min_distance = min_distance
        self.max_distance = max_distance
        self.max_order = max_order
        self.naming_pattern = naming_pattern
        self.exchange_type = exchange_type
        self.exchange_parameters = exchange_parameters
        self.direction_filter = direction_filter


    def exchanges(self, sites: list[Site], unit_cell: UnitCell):
        """Get the exchanges defined by this groups given the sites and symmetry provided.

        Parameters
        ----------
        sites : list of Site
            List of sites between which to generate exchanges.
        unit_cell : UnitCell
            Unit cell used to compute Cartesian bond vectors and distances.

        Returns
        -------
        list of Exchange
            Generated and filtered `Exchange` instances.
        """
        if self.bond > 0:
            # Try to improve this heuristic for the maximum distance
            distances = (0.0, np.sqrt(self.bond) * 2 * unit_cell.main_diagonal_length)
        else:
            distances = (self.min_distance, self.max_distance)

        abstract_exchanges = batch_exchanges(
            sites=sites,
            unit_cell=unit_cell,
            max_distance=distances[1],
            naming_pattern=default_naming_pattern if self.naming_pattern is None else self.naming_pattern,
            type_symbol=self.exchange_type.short_string)

        if self.bond > 0:
            dp = int(np.ceil(np.log10(1 / tolerances.BOND_TOL)))
            bond_distances = np.sort(np.unique(np.round([c.distance for c in abstract_exchanges], decimals=dp)))
            abstract_exchanges = [c for c in abstract_exchanges
                                  if np.abs(c.distance - bond_distances[self.bond-1]) < tolerances.BOND_TOL]

        kept_exchanges = []
        for exchange in abstract_exchanges:
            if self.max_order is not None and exchange.order > self.max_order:
                continue

            if distances[0] <= exchange.distance <= distances[1]:
                kept_exchanges.append(self.exchange_type(
                    name=exchange.name,
                    site_1=exchange.site_1,
                    site_2=exchange.site_2,
                    cell_offset=exchange.cell_offset,
                    **self.exchange_parameters))

        # Apply the filtering
        if self.direction_filter is not None:
            kept_exchanges = [exchange
                              for exchange in kept_exchanges
                              if self.direction_filter.accept(exchange.vector(unit_cell))]

        return kept_exchanges
