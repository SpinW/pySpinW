import numpy as np

from pyspinw.anisotropy import Anisotropy
from pyspinw.basis import site_rotations
from pyspinw.coupling import Coupling
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell

def alt_rotate(target_vector):
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
    def __init__(self, model_matrix: np.ndarray):
        self.model_matrix = model_matrix

class Selectable(Component):
    def __init__(self, render_id: int, model_matrix: np.ndarray):
        super().__init__(model_matrix)

        self.render_id = render_id
        self.is_selected = False
        self.is_hover = False
        self.is_group_selected = False
        self.is_group_hover = False


class RenderSite(Selectable):
    def __init__(self, render_id: int, site: LatticeSite, unit_cell: UnitCell):
        # Build model matrix

        rotation = site_rotations(np.array([site.base_moment]))
        translation = unit_cell.fractional_to_cartesian(site.ijk)

        model_matrix = np.zeros((4, 4), dtype=np.float32)
        model_matrix[3,3] = 1.0
        model_matrix[:3, :3] = -rotation
        model_matrix[:3, 3] = translation

        super().__init__(render_id, model_matrix)
        self.site = site


class RenderCoupling(Selectable):
    def __init__(self, render_id: int, coupling: Coupling, unit_cell: UnitCell):

        translation = unit_cell.fractional_to_cartesian(coupling.site_1.ijk)
        delta = unit_cell.fractional_to_cartesian(coupling.site_2.ijk + coupling.cell_offset.vector) - translation

        length = np.sqrt(np.sum(delta**2))
        print(length)

        rotation = alt_rotate(delta)
        print(delta)
        print(rotation)

        model_matrix = np.zeros((4, 4), dtype=np.float32)
        model_matrix[3,3] = 1.0
        model_matrix[:3, :3] = (rotation @ np.diag([1,1,length]) )
        model_matrix[:3, 3] = translation

        super().__init__(render_id, model_matrix)
        self.coupling = coupling

class RenderAnisotropy(Component):
    # Don't subclass selectable as we don't want to select anisotropies from the viewer

    def __init__(self, anisotropy: Anisotropy):

        super().__init__(np.eye(4))
        self.anisotropy = anisotropy

class RenderModel:
    def __init__(self, hamiltonian: Hamiltonian):

        self.hamiltonian = hamiltonian

        #
        # Build the data structure for rendering
        #

        self.expanded, site_mapping = self.hamiltonian.expand_with_mapping()

        render_id = 0

        self.render_map: dict[int, Component] = {}
        self.anisotropy_render_id_map: dict[int, int] = {} # We don't select anisotropies from the viewer, need this
        site_unique_id_to_render_id: dict[int, int] = {}

        self.sites = []
        for site in self.expanded.structure.sites:
            render_site = RenderSite(render_id, site, self.expanded.structure.unit_cell)
            self.sites.append(render_site)
            self.render_map[render_id] = render_site
            site_unique_id_to_render_id[site.unique_id] = render_id
            render_id += 1

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
            self.anisotropy_render_id_map[render_id] = site_unique_id_to_render_id[anisotropy.site.unique_id]
            render_id += 1

        #
        # Build the data structures for displaying
        #

