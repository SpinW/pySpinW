from collections import defaultdict

import numpy as np

from pyspinw.anisotropy import Anisotropy
from pyspinw.basis import site_rotations
from pyspinw.coupling import Coupling
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell

def alt_rotate(target_vector):
    """ Rotation matrix from (0,0,1) to the target vector direction """
    v = target_vector / np.sqrt(np.sum(target_vector**2))

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
    def __init__(self, render_id: int, site: LatticeSite, unit_cell: UnitCell, offset: tuple[int, int, int] | None):
        # Build model matrix

        rotation = site_rotations(np.array([site.base_moment]))
        translation = unit_cell.fractional_to_cartesian(site.ijk)

        model_matrix = np.zeros((4, 4), dtype=np.float32)
        model_matrix[3,3] = 1.0
        model_matrix[:3, :3] = -rotation
        model_matrix[:3, 3] = translation

        super().__init__(render_id, model_matrix)
        self.site = site
        self.offset = offset


class RenderCoupling(Selectable):
    """ Render information for each coupling """

    def __init__(self, render_id: int, coupling: Coupling, unit_cell: UnitCell):

        # The model this is designed for is a tube that goes from (0,0,0) to (0,0,1)

        translation = unit_cell.fractional_to_cartesian(coupling.site_1.ijk)
        delta = unit_cell.fractional_to_cartesian(coupling.site_2.ijk + coupling.cell_offset.vector) - translation

        length = np.sqrt(np.sum(delta**2))

        rotation = alt_rotate(delta)

        model_matrix = np.zeros((4, 4), dtype=np.float32)
        model_matrix[3,3] = 1.0
        model_matrix[:3, :3] = (rotation @ np.diag([1,1,length]) )
        model_matrix[:3, 3] = translation

        super().__init__(render_id, model_matrix)
        self.coupling = coupling

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
        # Build the data structures for displaying
        #
