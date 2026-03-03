import numpy as np

from pyspinw.gui.rendermodel import rotation_from_z
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.site import LatticeSite

class MinimisationConstraint:
    name = "<constraint base class>"

    def __repr__(self):
        return self.name

class FreeConstraint(MinimisationConstraint):
    name = "Free"

class FixedConstraint(MinimisationConstraint):
    name = "Fixed"

class Planar(MinimisationConstraint):
    name = "Planar"

    def __init__(self, axis=np.ndarray):
        self.axis = axis

    def __repr__(self):
        return f"{self.name}(axis={self.axis[0]},{self.axis[1]},{self.axis[2]})"

Free = FreeConstraint()
Fixed = FixedConstraint()

class ClassicalEnergyMinimisation:
    """
    Do a classical energy minimisation of the spin orientations

    See dev note 007_energy_minimisation.md
    """

    def __init__(self, hamiltonian: Hamiltonian, constraints: list[Constraint], field: np.ndarray):

        self.hamiltonian = Hamiltonian
        self.constraints = constraints
        self.field = field

        #
        # Gather the data we'll need for the calculation
        #

        # Split sites by constraint type

        sites = hamiltonian.structure.sites

        self.free_sites = []
        self.fixed_sites = []
        self.planar_sites = []
        self.planar_axes = []

        for i, (site, constraint) in enumerate(zip(sites, constraints)):
            match constraint:
                case FreeConstraint():
                    self.free_sites.append((i, site.unique_id))
                case FixedConstraint():
                    self.fixed_sites.append((i, site.unique_id))
                case PlanarConstraint(axis):
                    self.planar_sites.append((i, site.unique_id))
                    self.planar_axes.append(axis)

        self.n_free = len(self.free_sites)
        self.n_planar = len(self.planar_sites)
        self.n_fixed = len(self.fixed_sites)

        # Make lists of sites for each coupling, anisotropy
        self.site_to_coupling_side_1 = defaultdict(list)
        self.site_to_coupling_side_2 = defaultdict(list)
        self.site_to_anisotropy = defaultdict(list)

        for coupling in hamiltonian.couplings:
            self.site_to_coupling_side_1[coupling.site_1.unique_id].append(coupling)
            self.site_to_coupling_side_2[coupling.site_2.unique_id].append(coupling)

        for anisotropy in hamiltonian.anisotropies:
            self.site_to_anisotropy.append(anisotropy)

        # Magnetic field contributions
        self.field_contribution_vector = [] # Vectors v such that E_field = v.S = B g S
        for site in sites:
            self.field_contribution_vector.append(field @ site.g)

        # Get the derived quantities
        self.moments = np.array([site.base_moment for site in sites])
        self.magnitudes = np.sqrt(np.sum(self.moments**2, axis=1))
        self.directions = self.moments / self.magnitudes

        self._site_uid_to_index = {site.unique_id: i for i, site in enumerate(sites)}
        self.n_sites = len(sites)


    def iterate(self):

        # TODO: Modify account for supercells

        rotation_matrices = [rotation_from_z(direction) for direction in self.directions]


        forces_free_alpha = np.zeros((self.n_free,))
        forces_free_beta = np.zeros((self.n_free,))
        forces_planar = np.zeros((self.n_planar))

        for i, site_uid in self.free_sites:

            dS_dalpha = -rotation_matrices[:, 1] # m.(0, -1, 0)
            dS_dbeta = rotation_matrices[:, 0]   # m.(1,  0, 0)

            # Couplings
            for coupling in self.site_to_coupling_side_1[site_uid]:

                other_index = self._site_uid_to_index[coupling.site_2.unique_id]
                other_moment = self.moments[other_index, :]

                forces_free_alpha[i] += dS_dalpha @ coupling.coupling_matrix @ other_moment
                forces_free_beta[i] += dS_dbeta @ coupling.coupling_matrix @ other_moment

            for coupling in self.site_to_coupling_side_2[site_uid]:
                other_index = self._site_uid_to_index[coupling.site_1.unique_id]
                other_moment = self.moments[other_index, :]

                forces_free_alpha[i] += other_moment @ coupling.coupling_matrix @ dS_dalpha
                forces_free_beta[i] += other_moment @ coupling.coupling_matrix @ dS_dalpha

            # Anisotropies
            for anisotropy in self.site_to_anisotropy[site_uid]:
                pass

            # Field
            field_force_alpha = self.field_contribution_vector[i] @ dS_dalpha
            field_force_beta = self.field_contribution_vector[i] @ dS_dbeta
