""" Groups of couplings, a neater way of representing multiple, similar couplings

TODO: Currently broken - WIP

"""

import numpy as np
from numpy._typing import ArrayLike

from pyspinw.exchange import Exchange
from pyspinw.batch_couplings import batch_exchanges, default_naming_pattern
from pyspinw.symmetry.symmetry_settings import SymmetrySettings
from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.tolerances import tolerances


class AbstractCouplingGroup:
    """ Base class for coupling groups"""

    def couplings(self, sites: list[LatticeSite], symmetry: SymmetrySettings):
        """ Get the couplings defined by this groups given the sites and symmetry provided """
        raise NotImplementedError("couplings not available in base class")

class DirectionalityFilter:
    """ Base class for filtering couplings by direction"""

    def accept(self, vector: np.ndarray) -> bool:
        """ return True if in an allowed direction"""
        raise NotImplementedError("accept() not implemented in base class")

class InDirectionFilter(DirectionalityFilter):
    """ Selects vectors in a given direction """

    def __init__(self, direction: ArrayLike, max_dev_angle_deg: float = 0.01):
        direction = np.array(direction, dtype=float)
        self.direction = direction / np.sqrt(np.sum(direction**2))
        self.in_direction_dev_num = np.cos((np.pi/180) * max_dev_angle_deg)

    def accept(self, vector):
        """ DirectionalityFilter implementation: return True if in a similar direction"""
        sq_mag = np.sum(vector**2)
        if sq_mag == 0:
            return False
        else:
            return np.dot(self.direction, vector) / np.sqrt(sq_mag) > self.in_direction_dev_num


class BiDirectionFilter(DirectionalityFilter):
    """ Selects vectors in a given direction """

    def __init__(self, direction: ArrayLike, max_dev_angle_deg: float = 0.01):
        direction = np.array(direction, dtype=float)
        self.direction = direction / np.sqrt(np.sum(direction**2))
        self.in_direction_dev_num = np.cos((np.pi/180) * max_dev_angle_deg)

    def accept(self, vector):
        """ DirectionalityFilter implementation: return True if in a similar direction"""
        sq_mag = np.sum(vector**2)
        if sq_mag == 0:
            return False
        else:
            return np.abs(np.dot(self.direction, vector) / np.sqrt(sq_mag)) > self.in_direction_dev_num


class InPlaneFilter(DirectionalityFilter):
    """ Selects vectors in a given plane (specified by normal)"""

    def __init__(self, direction: ArrayLike, max_dev_angle_deg: float = 0.01):
        direction = np.array(direction)
        self.direction = direction / np.sqrt(np.sum(direction**2))
        self.in_plane_dev_num = np.sin((np.pi/180) * max_dev_angle_deg)

    def accept(self, vector):
        """ DirectionalityFilter implementation: return True if in the plane normal to provided vector"""
        sq_mag = np.sum(vector ** 2)
        if sq_mag == 0:
            return False
        else:
            return np.abs(np.dot(self.direction, vector) / np.sqrt(sq_mag)) < self.in_plane_dev_num


class ExchangeGroup:
    """ Class representing the batch creation of couplings"""

    def __init__(self,
                 name: str,
                 bond: int,
                 min_distance: float,
                 max_distance: float,
                 max_order: int | None,
                 naming_pattern: str | None,
                 exchange_type: type[Exchange],
                 coupling_parameters: dict,
                 direction_filter: DirectionalityFilter | None):

        self.name = name
        self.bond = bond
        self.min_distance = min_distance
        self.max_distance = max_distance
        self.max_order = max_order
        self.naming_pattern = naming_pattern
        self.exchange_type = exchange_type
        self.coupling_parameters = coupling_parameters
        self.direction_filter = direction_filter


    def exchanges(self, sites: list[LatticeSite], unit_cell: UnitCell):
        """ Get the couplings defined by this groups given the sites and symmetry provided """
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
                    **self.coupling_parameters))

        # Apply the filtering
        if self.direction_filter is not None:
            kept_exchanges = [exchange
                              for exchange in kept_exchanges
                              if self.direction_filter.accept(exchange.vector(unit_cell))]

        return kept_exchanges


