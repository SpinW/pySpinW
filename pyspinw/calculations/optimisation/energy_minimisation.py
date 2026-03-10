from collections import defaultdict

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

alpha_m = np.array([0,-1,0], dtype=float)
beta_m = np.array([1,0,0], dtype=float)

class ClassicalEnergyMinimisation:
    """
    Do a classical energy minimisation of the spin orientations

    See dev note 007_energy_minimisation.md
    """

    def __init__(self, hamiltonian: Hamiltonian, constraints: list[MinimisationConstraint], field: np.ndarray):

        self.hamiltonian = hamiltonian
        self.constraints = constraints
        self.field = field

        #
        # Gather the data we'll need for the calculation
        #

        # Split sites by constraint type

        sites = hamiltonian.structure.sites
        self.n_sites = len(sites)

        self.free_sites = []
        self.fixed_sites = []
        self.planar_sites = []
        self.planar_axes = []

        self.is_free = np.zeros((self.n_sites,), dtype=bool)
        self.is_fixed = np.zeros((self.n_sites,), dtype=bool)
        self.is_planar = np.zeros((self.n_sites,), dtype=bool)

        for i, (site, constraint) in enumerate(zip(sites, constraints)):
            match constraint:

                case FreeConstraint():
                    self.free_sites.append((i, site.unique_id))
                    self.is_free[i] = True

                case FixedConstraint():
                    self.fixed_sites.append((i, site.unique_id))
                    self.is_fixed[i] = True

                case Planar(axis):
                    self.planar_sites.append((i, site.unique_id))
                    self.planar_axes.append(axis)
                    self.is_planar[i] = True


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
            self.site_to_anisotropy[anisotropy.site.unique_id].append(anisotropy)

        # Magnetic field contributions
        self.field_contribution_vector = [] # Vectors v such that E_field = v.S = B g S
        for site in sites:
            self.field_contribution_vector.append(field @ site.g)

        # Get the derived quantities
        self.moments = np.array([site.base_moment for site in sites])
        self.magnitudes = np.sqrt(np.sum(self.moments**2, axis=1))

        self._site_uid_to_index = {site.unique_id: i for i, site in enumerate(sites)}

    def jitter(self, max_angular_change_per_direction):
        pass

    def energy(self):
        energy = 0.0
        for coupling in self.hamiltonian.couplings:
            site_1_moment = self.moments[self._site_uid_to_index[coupling.site_1.unique_id], :]
            site_2_moment = self.moments[self._site_uid_to_index[coupling.site_2.unique_id], :]

            energy += site_1_moment @ coupling.coupling_matrix @ site_2_moment

        # TODO: Anisotropies and fields

        return energy

    def iterate(self, step_size_factor=0.01):

        # TODO: Modify account for supercells

        rotation_matrices = [rotation_from_z(moment) for moment in self.moments]


        forces_free_alpha = np.zeros((self.n_free,))
        forces_free_beta = np.zeros((self.n_free,))
        forces_planar = np.zeros((self.n_planar))

        for param_index, (site_index, site_uid) in enumerate(self.free_sites):

            # dS_dalpha = -rotation_matrices[i][:, 1] # m.(0, -1, 0)
            # dS_dbeta = rotation_matrices[i][:, 0]   # m.(1,  0, 0)

            dS_dalpha = rotation_matrices[site_index] @ alpha_m # m.(0, -1, 0)
            dS_dbeta = rotation_matrices[site_index] @ beta_m   # m.(1,  0, 0)


            # Couplings
            for coupling in self.site_to_coupling_side_1[site_uid]:

                other_index = self._site_uid_to_index[coupling.site_2.unique_id]
                other_moment = self.moments[other_index, :]

                forces_free_alpha[site_index] -= dS_dalpha @ coupling.coupling_matrix @ other_moment
                forces_free_beta[site_index] -= dS_dbeta @ coupling.coupling_matrix @ other_moment

            for coupling in self.site_to_coupling_side_2[site_uid]:
                other_index = self._site_uid_to_index[coupling.site_1.unique_id]
                other_moment = self.moments[other_index, :]

                forces_free_alpha[param_index] -= other_moment @ coupling.coupling_matrix @ dS_dalpha
                forces_free_beta[param_index] -= other_moment @ coupling.coupling_matrix @ dS_dbeta

            # Anisotropies
            for anisotropy in self.site_to_anisotropy[site_uid]:
                pass

            # Field
            field_force_alpha = self.field_contribution_vector[site_index] @ dS_dalpha
            field_force_beta = self.field_contribution_vector[site_index] @ dS_dbeta

            forces_free_alpha[param_index] -= field_force_alpha
            forces_free_beta[param_index] -= field_force_beta

        #
        # print("alpha forces:", forces_free_alpha)
        # print("beta forces:", forces_free_beta)

        # Move in direction of force
        # If you want a physical interpretation, this is critically damped movement, where step_size_factor is dt
        # AKA 90s game engine physics, AKA Aristotlean physics, F = mv
        #
        # As we're in a coordinate system around current location, alpha = delta alpha

        alpha = step_size_factor * np.array(forces_free_alpha)
        beta = step_size_factor * np.array(forces_free_beta)

        # print("alpha change:", alpha)
        # print("beta change:", beta)

        # Get the moments before rotation

        new_moments = np.zeros((self.n_sites, 3), dtype=float)
        new_moments[:, 2] = 1.0

        cos_beta = np.cos(beta)

        new_moments[self.is_free, :] = np.array([
            np.sin(beta),
            -np.sin(alpha) * cos_beta,
            np.cos(alpha) * cos_beta
        ]).T

        for site_index, _ in self.free_sites:
            self.moments[site_index] = rotation_matrices[site_index] @ new_moments[site_index, :]