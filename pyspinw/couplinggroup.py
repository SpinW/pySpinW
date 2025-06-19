import numpy as np

from pyspinw._base import Coupling
from pyspinw.batch_couplings import batch_couplings
from pyspinw.gui.symmetry_settings import SymmetrySettings
from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell


class AbstractCouplingGroup:


    def couplings(self, sites: list[LatticeSite], symmetry: SymmetrySettings):
        raise NotImplementedError("couplings not available in base class")

class DirectionalityFilter:
    def accept(self, vector: np.ndarray):
        raise NotImplementedError("__call__ not implemented in base class")

class InDirectionFilter(DirectionalityFilter):
    def __init__(self, direction: np.ndarray, max_dev_angle_deg: float):
        self.direction = direction / np.sqrt(np.sum(direction**2))
        self.in_direction_dev_num = np.cos((np.pi/180) * max_dev_angle_deg)

    def accept(self, vector):
        sq_mag = np.sum(vector**2)
        if sq_mag == 0:
            return False
        else:
            return np.dot(self.direction, vector) / np.sqrt(sq_mag) > self.in_direction_dev_num

class InPlaneFilter(DirectionalityFilter):
    def __init__(self, direction: np.ndarray, max_dev_angle_deg: float):
        self.direction = direction / np.sqrt(np.sum(direction**2))
        self.in_plane_dev_num = np.sin((np.pi/180) * max_dev_angle_deg)

    def accept(self, vector):
        sq_mag = np.sum(vector ** 2)
        if sq_mag == 0:
            return False
        else:
            return np.dot(self.direction, vector) / np.sqrt(sq_mag) < self.in_plane_dev_num


class CouplingGroup():
    """ Class representing the batch creation of couplings"""

    def __init__(self,
                 name: str,
                 site_list_indices: list[int],
                 min_distance: float,
                 max_distance: float,
                 naming_pattern: str,
                 coupling_type: type[Coupling],
                 parameters: dict,
                 direction_filter: DirectionalityFilter):

        self.name = name
        self.min_distance = min_distance
        self.max_distance = max_distance
        self.naming_pattern = naming_pattern
        self.coupling_type = coupling_type
        self.parameters = parameters
        self.direction_filter = direction_filter


    def couplings(self, sites: list[LatticeSite], unit_cell: UnitCell):

        abstract_couplings = batch_couplings(
            sites=self.sites,
            unit_cell=self.unit_cell,
            max_distance=self.max_distance,
            naming_pattern=self.format_string,
            type_symbol=self.short_string)

        parameter_dict = {parameter: self.field_widget.get_value(parameter)
                          for parameter in coupling_type.parameters}

        kept_couplings = []
        for coupling in abstract_couplings:
            if parameters.max_order is not None and coupling.order > parameters.max_order:
                continue

            if parameters.min_distance <= coupling.distance <= parameters.max_distance:
                kept_couplings.append(coupling_type(
                    name=coupling.name,
                    site_1=coupling.site_1,
                    site_2=coupling.site_2,
                    cell_offset=coupling.cell_offset,
                    **parameter_dict))

        return kept_couplings

class SetCouplingTypeAndParameters(AbstractCouplingGroup):
    def __init__(self, parent: AbstractCouplingGroup, indices: list[int] | None, coupling_type: str | None = None, parameters: dict = {}):
        self.parent = parent

class SetName(AbstractCouplingGroup):
    def __init__(self, parent: AbstractCouplingGroup, *index_name_pairs: tuple[int, str]):
        self.parent = parent

class Remove(AbstractCouplingGroup):
    def __init__(self, parent: AbstractCouplingGroup, *indices: int):
        self.parent = parent

        