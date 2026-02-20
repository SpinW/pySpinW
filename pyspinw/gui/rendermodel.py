from collections import defaultdict

import numpy as np

from pyspinw.anisotropy import Anisotropy
from pyspinw.basis import site_rotations
from pyspinw.coupling import Coupling
from pyspinw.gui.edge_cases import add_extra_edge_lines, add_extra_edge_points
from pyspinw.gui.wrap_line import split_and_wrap_line_segment
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell


def rotation_from_z(target_vector):
    """ Rotation matrix from (0,0,1) to the target vector direction """
    mag_sq = np.sum(target_vector**2)

    if mag_sq < 1e-9:
        return np.eye(3)

    v = target_vector / np.sqrt(mag_sq)

    x,y,z = v

    if z < 1e-9 - 1:
        return np.array([[-1,0,0],[0,1,0],[0,0,-1]])

    return np.array([
        [1-x**2 / (1+z), -x*y / (1+z), x],
        [-x*y / (1+z), 1 - y**2 / (1+z), y],
        [-x, -y, z]
    ])


class Component:
    """ Base class for components """
    def __init__(self, model_matrix: np.ndarray):
        self.model_matrix = model_matrix

class Selectable(Component):
    """ Components that can be selected in the viewer """
    def __init__(self, render_id: int, model_matrix: np.ndarray):
        super().__init__(model_matrix)

        self.render_id = render_id


class RenderSite(Selectable):

    """ Render information for each site """

    @staticmethod
    def site_model_matrix(position, moment, unit_cell):

        rotation = rotation_from_z(moment)
        translation = unit_cell.fractional_to_cartesian(position)

        model_matrix = np.zeros((4, 4), dtype=np.float32)
        model_matrix[3,3] = 1.0
        model_matrix[:3, :3] = rotation
        model_matrix[:3, 3] = translation

        return model_matrix

    def __init__(self, render_id: int, site: LatticeSite, unit_cell: UnitCell, offset: tuple[int, int, int] | None):
        # Build model matrix

        position, moment = site.ijk, site.base_moment
        model_matrix = self.site_model_matrix(position, moment, unit_cell)

        super().__init__(render_id, model_matrix)
        self.site = site
        self.offset = offset

        self.pretty_render_model_matrices = [self.site_model_matrix(new_position, moment, unit_cell)
                                             for new_position in add_extra_edge_points(position)]


    def model_matrices(self, pretty: bool) -> list[np.ndarray]:

        if pretty:
            return self.pretty_render_model_matrices
        else:
            return [self.model_matrix]


class RenderCoupling(Selectable):
    """ Render information for each coupling """

    @staticmethod
    def segment_model_matrix(a, b, unit_cell):

        # The model this is designed for is a tube that goes from (0,0,0) to (0,0,1)

        translation = unit_cell.fractional_to_cartesian(a)
        delta = unit_cell.fractional_to_cartesian(b) - translation

        length = np.sqrt(np.sum(delta**2))

        rotation = rotation_from_z(delta)

        model_matrix = np.zeros((4, 4), dtype=np.float32)
        model_matrix[3,3] = 1.0
        model_matrix[:3, :3] = (rotation @ np.diag([1,1,length]) )
        model_matrix[:3, 3] = translation

        return model_matrix

    def __init__(self, render_id: int, coupling: Coupling, unit_cell: UnitCell):

        self.coupling = coupling

        p1 = coupling.site_1.ijk
        p2 = coupling.site_2.ijk + coupling.cell_offset.vector

        model_matrix = self.segment_model_matrix(p1, p2, unit_cell)

        super().__init__(render_id, model_matrix)

        # Create prettified model matrices
        _, sections = split_and_wrap_line_segment(p1, p2)
        sections = add_extra_edge_lines(sections)
        self.split_model_matrices = [self.segment_model_matrix(s1, s2, unit_cell) for s1, s2 in sections]

    def model_matrices(self, pretty: bool) -> list[np.ndarray]:
        """Either the default matrix in a list, or a list of the prettified matrices"""

        if pretty:
            return self.split_model_matrices
        else:
            return [self.model_matrix]


class RenderAnisotropy(Component):
    """ Render information for anisotropies"""

    # Don't subclass selectable as we don't want to select anisotropies from the viewer

    def __init__(self, anisotropy: Anisotropy):

        super().__init__(np.eye(4))
        self.anisotropy = anisotropy

class RenderModel:
    """ Model of the hamiltonian, contains lots of derived information for rendering and text representation"""

    def __init__(self, hamiltonian: Hamiltonian):

        self.hamiltonian = hamiltonian

        self.expanded, site_mapping, self.coupling_mapping, self.anisotropy_mapping = \
            self.hamiltonian._expand_with_mapping()

        self.site_expanded_uid_to_original_uid: dict[int, int] = {}
        self.site_original_uid_to_expanded_uid: defaultdict[int, list[int]] = defaultdict(list)

        self.original_sites = self.hamiltonian.structure.sites

        for (parent_unique_id, offset), child_site in site_mapping.items():
            self.site_original_uid_to_expanded_uid[parent_unique_id].append(child_site.unique_id)
            self.site_expanded_uid_to_original_uid[child_site.unique_id] = parent_unique_id

        self.site_expanded_uid_to_original_index_and_offset = \
            {site.unique_id: data for data, site in site_mapping.items()}

        #
        # Build the data structure for rendering
        #

        render_id = 1

        self.render_map: dict[int, Component] = {}
        self.expanded_site_unique_id_to_render_id: dict[int, int] = {}

        self.original_site_unique_id_to_index: dict[int, int] = {}
        for index, site in enumerate(self.hamiltonian.structure.sites):
            self.original_site_unique_id_to_index[site.unique_id] = index

        self.expanded_index_to_original_index = []

        #
        # Create render data for sites
        #
        self.sites: list[RenderSite] = []
        for site in self.expanded.structure.sites:

            parent_index, offset = self.site_expanded_uid_to_original_index_and_offset[site.unique_id]
            render_site = RenderSite(render_id, site, self.expanded.structure.unit_cell, offset)

            self.sites.append(render_site)
            self.render_map[render_id] = render_site
            self.expanded_site_unique_id_to_render_id[site.unique_id] = render_id

            parent_unique_id = self.site_expanded_uid_to_original_uid[site.unique_id]

            self.expanded_index_to_original_index.append(self.original_site_unique_id_to_index[parent_unique_id])

            render_id += 1

        self.original_index_to_expanded = defaultdict(list)
        for expanded, original in enumerate(self.expanded_index_to_original_index):
            self.original_index_to_expanded[original].append(expanded)

        #
        # Create render data for couplings
        #

        self.couplings = []
        for coupling in self.expanded.couplings:
            render_coupling = RenderCoupling(render_id, coupling, self.expanded.structure.unit_cell)
            self.couplings.append(render_coupling)
            self.render_map[render_id] = render_coupling
            render_id += 1

        self.anisotropies = []
        for anisotropy in self.expanded.anisotropies:
            render_anisotropy = RenderAnisotropy(anisotropy)
            self.anisotropies.append(render_anisotropy)
            self.render_map[render_id] = render_anisotropy
            render_id += 1

        #
        # Unit cell axes
        #

        translation_matrix = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,.5],[0,0,0,1]], dtype=np.float32)

        # get orthogonal axes transforms

        unit_cell_orthobasis = hamiltonian.structure.unit_cell._xyz_moments

        self.unit_cell_orthobasis_transforms = []
        for direction in unit_cell_orthobasis:
            rotation_matrix = np.eye(4, dtype=np.float32)
            rotation_matrix[:3, :3] = rotation_from_z(direction)

            self.unit_cell_orthobasis_transforms.append(rotation_matrix @ translation_matrix)

        # get normal xyz

        unit_cell_basis = hamiltonian.structure.unit_cell._xyz

        self.unit_cell_axes_transforms = []
        for direction in unit_cell_basis:
            rotation_matrix = np.eye(4, dtype=np.float32)
            rotation_matrix[:3, :3] = rotation_from_z(direction)

            self.unit_cell_axes_transforms.append(rotation_matrix @ translation_matrix)
