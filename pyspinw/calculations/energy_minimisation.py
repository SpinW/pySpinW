""" Classical minimisation of spin orientations """

from collections import defaultdict
from enum import Enum
from abc import ABC, abstractmethod

import numpy as np
from numpy._typing import ArrayLike

from pyspinw.symmetry.supercell import CommensurateSupercell
from pyspinw.util import triple_product_matrix, rotation_matrix, rotation_from_z



class MinimisationConstraintGenerator(ABC):
    """ Base class for constraint specifications that can be applied to all sites"""

    @abstractmethod
    def generate(self, n_constraints: int):
        """ Generate full list of constraints"""

class MinimisationConstraint:
    """ Base class for constraint specifications applied to a single site"""

    name = "<constraint base class>"

    def __repr__(self):
        return self.name

class FreeConstraint(MinimisationConstraint, MinimisationConstraintGenerator):
    """ Specify that moment should be unconstrained """

    name = "Free"

    def generate(self, n_constraints: int):
        """ Generate a list of constraints """
        return [self for _ in range(n_constraints)]


class FixedConstraint(MinimisationConstraint):
    """ Specify that moment should be fixed """

    name = "Fixed"

class Planar(MinimisationConstraint, MinimisationConstraintGenerator):
    """ Specify that moment should be constrained to plane """

    name = "Planar"

    __match_args__ = ('axis', )

    def __init__(self, axis: ArrayLike):
        self.axis = np.array(axis, dtype=float)
        self.axis /= np.sqrt(np.sum(self.axis**2))

    def generate(self, n_constraints: int):
        """ Generate a list of constraints """
        return [self for _ in range(n_constraints)]

    def __repr__(self):
        return f"{self.name}(axis={self.axis[0]},{self.axis[1]},{self.axis[2]})"

Free = FreeConstraint()
Fixed = FixedConstraint()

class OneFixed(MinimisationConstraintGenerator):
    """ Constraint Generator"""

    def __init__(self, index:int=0):
        self.index = index

    def generate(self, n_constraints: int):
        """ Generate a list of constraints """
        return [Fixed if self.index == i else Free for i in range(n_constraints)]

class InitialRandomisation(Enum):
    """ Enum for setting the method of initial randomisation """

    NONE = "none"
    JITTER = "jitter"
    RANDOMISED = "randomised"


alpha_m = np.array([0,-1,0], dtype=float)
beta_m = np.array([1,0,0], dtype=float)

class ClassicalEnergyMinimisation:
    """Do a classical energy minimisation of the spin orientations

    See dev note 007_energy_minimisation.md
    """

    def __init__(self,
                 hamiltonian: "Hamiltonian",
                 constraints: list[MinimisationConstraint] | MinimisationConstraintGenerator = Free,
                 field: ArrayLike | None = None,
                 seed: int | None = None):

        # Check that the supercell is commensurate

        if not isinstance(hamiltonian.structure.supercell, CommensurateSupercell):
            raise TypeError("Classical energy minimisation only implemented for commensurate structures")

        self.supercell: CommensurateSupercell = hamiltonian.structure.supercell

        # Check / normalise parameters

        n_sites = len(hamiltonian.structure.sites)

        if isinstance(constraints, list):
            if len(constraints) != n_sites:
                raise ValueError("If constraints are a list, there should be one for each site")

        elif isinstance(constraints, MinimisationConstraintGenerator):
            constraints = constraints.generate(n_sites)

        else:
            raise TypeError("Expected constraints to be a MinimisationConstraintGenerator or "
                            "list of MinimisationConstraints")


        if field is None:
            field = np.zeros((3, ))
        else:
            field = np.array(field)

        if field.shape != (3, ):
            raise ValueError("Expected field to be a length 3 vector")

        # Set up basic fields

        self.hamiltonian = hamiltonian
        self.constraints = constraints
        self.field = field

        self.rng = np.random.default_rng(seed)

        #
        # Gather the data we'll need for the calculation
        #

        # Split sites by constraint type

        self.sites = hamiltonian.structure.sites
        self.n_sites = len(self.sites)
        self.n_components = hamiltonian.structure.supercell.n_components()

        self.free_sites = []
        self.fixed_sites = []
        self.planar_sites = []
        self.planar_axes = []

        self.is_free = np.zeros((self.n_sites,), dtype=bool)
        self.is_fixed = np.zeros((self.n_sites,), dtype=bool)
        self.is_planar = np.zeros((self.n_sites,), dtype=bool)

        for site_index, (site, constraint) in enumerate(zip(self.sites, constraints)):
            match constraint:

                case FreeConstraint():
                    self.free_sites.append((site_index, site.unique_id))
                    self.is_free[site_index] = True

                case FixedConstraint():
                    self.fixed_sites.append((site_index, site.unique_id))
                    self.is_fixed[site_index] = True

                case Planar(axis):
                    self.planar_sites.append((site_index, site.unique_id))
                    self.planar_axes.append(axis)
                    self.is_planar[site_index] = True


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
        for site in self.sites:
            self.field_contribution_vector.append(field @ site.g)

        # Get the derived quantities
        self.moment_data = np.array([site.moment_data for site in self.sites])

        assert self.moment_data.shape == (self.n_sites, self.n_components, 3)

        self.magnitudes = np.sqrt(np.sum(self.moment_data ** 2, axis=2)).reshape(-1, self.n_components, 1)

        assert self.magnitudes.shape == (self.n_sites, self.n_components, 1)

        self._site_uid_to_index = {site.unique_id: i for i, site in enumerate(self.sites)}

    @staticmethod
    def _random_orientations(rng, shape, dim=3):
        """ Make random gaussian vectors, avoiding 0,0,0..."""
        n = 1
        for m in shape:
            n *= m

        output = rng.normal(0, 1, (n, dim))

        # Repeatedly replace any zero rows
        zero_rows = np.all(output == 0, axis=1)
        n_bad = np.sum(np.array(zero_rows, int))

        while n_bad > 0:
            output[zero_rows, :] = rng.normal(0, 1, (n_bad, dim))
            zero_rows = np.all(output == 0, axis=1)
            n_bad = np.sum(np.array(zero_rows, int))

        output /= np.sqrt(np.sum(output**2, axis=1)).reshape(-1, 1)

        return output.reshape(*(shape + (dim, ))) # reshape to expected shape

    def jitter(self, jitter_size_rad=0.1):
        """ Jitter method applies a movement of fixed size (radians) in a random direction """
        #
        # Free rotations
        #

        rotation_matrices = [[rotation_from_z(self.moment_data[i, j, :])
                              for i, uid in self.free_sites]
                              for j in range(self.n_components)]

        # random direction in angle
        angles = (2*np.pi) * self.rng.random((self.n_free, ))
        cos_angles = np.cos(angles)
        sin_angles = np.sin(angles)

        cos_jitter = np.cos(jitter_size_rad)
        sin_jitter = np.sin(jitter_size_rad)

        unrotated_moments = np.empty((self.n_free, self.n_components, 3), dtype=float)
        for component_index in range(self.n_components):
            unrotated_moments[:, component_index, 0] = sin_jitter * cos_angles
            unrotated_moments[:, component_index, 1] = sin_jitter * sin_angles
            unrotated_moments[:, component_index, 2] = cos_jitter

        unrotated_moments *= self.magnitudes[self.is_free, :]

        for param_index, (site_index, _) in enumerate(self.free_sites):
            for component_index in range(self.n_components):
                self.moment_data[site_index, component_index, :] = \
                    rotation_matrices[component_index][param_index] @ unrotated_moments[param_index, component_index, :]

        #
        # Planar moments
        #

        for (site_index, _), axis in zip(self.planar_sites, self.planar_axes):
            for component_index in range(self.n_components):
                angle = jitter_size_rad if self.rng.random() > 0.5 else -jitter_size_rad
                self.moment_data[site_index, component_index, :] = \
                    rotation_matrix(angle, axis) @ self.moment_data[site_index, component_index, :]

    def randomise(self):
        """ Randomise method chooses random directions for spins"""
        # Do Free sites using Gaussian method

        random_orientations = self._random_orientations(self.rng, (self.n_free, self.n_components), 3)

        self.moment_data[self.is_free, :, :] = random_orientations
        self.moment_data[self.is_free, :, :] *= self.magnitudes[self.is_free, :, :]



        # Do planar sites by choosing a random number in [0,2pi)

        for param_index, ((site_index, site_uid), axis) in enumerate(zip(self.planar_sites, self.planar_axes)):
            for component_index in range(self.n_components):

                random_angle = (2 * np.pi) * self.rng.random()

                random_z_perp = np.array([
                    np.sin(random_angle),
                    np.cos(random_angle),
                    0.0]) * self.magnitudes[site_index, component_index, 0]

                self.moment_data[site_index, component_index, :] = rotation_from_z(axis) @ random_z_perp


    def energy(self):
        """ Energy of the current moments according to the hamiltonian, per unit cell """
        supercell = self.hamiltonian.structure.supercell

        energy = 0.0

        for cell in supercell.cells():

            # Exchanges
            for coupling in self.hamiltonian.couplings:
                site_1_moment_data = self.moment_data[self._site_uid_to_index[coupling.site_1.unique_id], :, :]
                site_2_moment_data = self.moment_data[self._site_uid_to_index[coupling.site_2.unique_id], :, :]

                site_1_moment = supercell.moment_calculation(site_1_moment_data, cell)
                site_2_moment = supercell.moment_calculation(
                                    site_2_moment_data,
                                    supercell.wrap_sum(cell, coupling.cell_offset))

                energy += site_1_moment @ coupling.coupling_matrix @ site_2_moment

            # Anisotropies
            for anisotropy in self.hamiltonian.anisotropies:
                moment_data = self.moment_data[self._site_uid_to_index[anisotropy.site.unique_id], :, :]

                moment = supercell.moment_calculation(moment_data, cell)

                energy += moment @ anisotropy.anisotropy_matrix @ moment

            # Field contribution
            for site_index in range(self.n_sites):
                moment_data = self.moment_data[site_index, :, :]

                moment = supercell.moment_calculation(moment_data, cell)

                energy += np.dot(self.field_contribution_vector[site_index], moment)

        return energy / supercell.n_cells()

    def minimise(self, rtol=1e-10, atol=1e-12, max_iters=1000,
                 initial_randomisation: InitialRandomisation | str = InitialRandomisation.JITTER,
                 verbose=False):
        """ Automatically do the minimisation and stop based on energy convergence """
        if isinstance(initial_randomisation, str):
            initial_randomisation = InitialRandomisation(initial_randomisation.lower())

        match initial_randomisation:
            case InitialRandomisation.NONE:
                pass
            case InitialRandomisation.JITTER:
                if verbose:
                    print("Jittering")
                self.jitter()
            case InitialRandomisation.RANDOMISED:
                if verbose:
                    print("Randomising")
                self.randomise()

        #
        # One iteration before main loop to get the initial energy change
        #

        last_energy = self.energy()
        self.iterate()
        this_energy = self.energy()

        start_delta_energy = this_energy - last_energy

        #
        # Main loop
        #

        for i in range(max_iters-1):
            last_energy = this_energy
            self.iterate()
            this_energy = self.energy()

            delta_energy = this_energy - last_energy

            if np.abs(delta_energy) < rtol * np.abs(start_delta_energy):
                if verbose:
                    print(f"Converged to E={this_energy} after {i+1} iterations "
                          f"(dE = {delta_energy} < {rtol} x {start_delta_energy})")

                break

            if np.abs(delta_energy) < atol:
                if verbose:
                    print(
                        f"Converged to E={this_energy} after {i + 1} iterations (dE = {delta_energy} < {atol})")

                break


        else:
            if verbose:
                print(f"Failed to converge after {max_iters} iterations")

        # TODO: Return new Hamiltonian



    def iterate(self, step_size_factor=0.1):
        """ Perform one step of gradient descent """
        # Free sites

        rotation_matrices = [[rotation_from_z(self.moment_data[site_index, component_index, :])
                                for component_index in range(self.n_components)]
                              for site_index, _ in self.free_sites]

        forces_free_alpha = np.zeros((self.n_free, self.n_components))
        forces_free_beta = np.zeros((self.n_free, self.n_components))

        for param_index, (site_index, site_uid) in enumerate(self.free_sites):
            for component_index in range(self.n_components):

                dT_dalpha = rotation_matrices[param_index][component_index] @ alpha_m  # m.(0, -1, 0)
                dT_dbeta = rotation_matrices[param_index][component_index] @ beta_m  # m.(1,  0, 0)

                for cell in self.supercell.cells():


                    # Couplings

                    for coupling in self.site_to_coupling_side_2[site_uid]:
                        # This side of the couplings has the offset applied

                        offset_cell = self.supercell.wrap_sum(cell, coupling.cell_offset)

                        df_dT = self.supercell.moment_derivative(component_index, offset_cell)

                        dS_dalpha = df_dT @ dT_dalpha
                        dS_dbeta = df_dT @ dT_dbeta

                        # Get other moment in supercell
                        other_index = self._site_uid_to_index[coupling.site_1.unique_id]
                        other_moment_data = self.moment_data[other_index, :, :]
                        other_moment = self.supercell.moment_calculation(other_moment_data, cell)

                        # Get forces
                        forces_free_alpha[param_index, component_index] -= \
                            other_moment @ coupling.coupling_matrix @ dS_dalpha

                        forces_free_beta[param_index, component_index] -= \
                            other_moment @ coupling.coupling_matrix @ dS_dbeta

                    # Calculate the derivative of the spin with respect to the parameters,
                    #  this depended on the cell for the other couplings, but is the same for
                    #  everything going forward
                    df_dT = self.supercell.moment_derivative(component_index, cell)

                    dS_dalpha = df_dT @ dT_dalpha
                    dS_dbeta = df_dT @ dT_dbeta


                    for coupling in self.site_to_coupling_side_1[site_uid]:

                        # The side without the offset applied

                        # Get other moment in supercell
                        other_index = self._site_uid_to_index[coupling.site_2.unique_id]
                        other_moment_data = self.moment_data[other_index, :, :]
                        other_cell = self.supercell.wrap_sum(cell, coupling.cell_offset)
                        other_moment = self.supercell.moment_calculation(other_moment_data, other_cell)

                        # Get the forces

                        forces_free_alpha[site_index, component_index] -= \
                            dS_dalpha @ coupling.coupling_matrix @ other_moment

                        forces_free_beta[site_index, component_index] -= \
                            dS_dbeta @ coupling.coupling_matrix @ other_moment

                    # Anisotropies
                    for anisotropy in self.site_to_anisotropy[site_uid]:
                        # dE = m.A.dm + (dm.A.m).T (.T not needed when using 1D arrays)
                        current_moment_data = self.moment_data[param_index, :, :]
                        current_moment = self.supercell.moment_calculation(current_moment_data, cell)

                        forces_free_alpha[param_index, component_index] -= \
                            current_moment @ anisotropy.anisotropy_matrix @ dS_dalpha

                        forces_free_alpha[param_index, component_index] -= \
                            dS_dalpha @ anisotropy.anisotropy_matrix @ current_moment

                        forces_free_beta[param_index, component_index] -= \
                            current_moment @ anisotropy.anisotropy_matrix @ dS_dbeta

                        forces_free_beta[param_index, component_index] -= \
                            dS_dbeta @ anisotropy.anisotropy_matrix @ current_moment

                    # Field
                    field_force_alpha = self.field_contribution_vector[site_index] @ dS_dalpha
                    field_force_beta = self.field_contribution_vector[site_index] @ dS_dbeta

                    forces_free_alpha[param_index, component_index] -= field_force_alpha
                    forces_free_beta[param_index, component_index] -= field_force_beta

        # Planar sites

        forces_planar = np.zeros((self.n_planar, self.n_components))

        for param_index, ((site_index, site_uid), axis) in enumerate(zip(self.planar_sites, self.planar_axes)):
            for component_index in range(self.n_components):

                dT_dtheta = triple_product_matrix(-axis) @ self.moment_data[site_index, component_index, :]

                for cell in self.supercell.cells():

                    for coupling in self.site_to_coupling_side_2[site_uid]:

                        offset_cell = self.supercell.wrap_sum(cell, coupling.cell_offset)

                        dS_dT = self.supercell.moment_derivative(component_index, offset_cell)
                        dS_dtheta = dS_dT @ dT_dtheta

                        # Get other moment in supercell
                        other_index = self._site_uid_to_index[coupling.site_1.unique_id]
                        other_moment_data = self.moment_data[other_index, :, :]
                        other_moment = self.supercell.moment_calculation(other_moment_data, cell)

                        forces_planar[param_index] -= other_moment @ coupling.coupling_matrix @ dS_dtheta

                    dS_dT = self.supercell.moment_derivative(component_index, cell)
                    dS_dtheta = dS_dT @ dT_dtheta

                    # Couplings
                    for coupling in self.site_to_coupling_side_1[site_uid]:

                        other_index = self._site_uid_to_index[coupling.site_2.unique_id]
                        other_moment_data = self.moment_data[other_index, :, :]
                        other_cell = self.supercell.wrap_sum(cell, coupling.cell_offset)
                        other_moment = self.supercell.moment_calculation(other_moment_data, other_cell)

                        forces_planar[site_index] -= dS_dtheta @ coupling.coupling_matrix @ other_moment

                    # Anisotropies
                    for anisotropy in self.site_to_anisotropy[site_uid]:
                        # dE = m.A.dm + (dm.A.m).T (.T not needed when using 1D arrays)

                        current_moment_data = self.moment_data[param_index, :, :]
                        current_moment = self.supercell.moment_calculation(current_moment_data, cell)

                        forces_planar[param_index, component_index] -= \
                            current_moment @ anisotropy.anisotropy_matrix @ dS_dtheta

                        forces_planar[param_index, component_index] -= \
                            dS_dtheta @ anisotropy.anisotropy_matrix @ current_moment


                    # Field
                    forces_planar[param_index, component_index] -= \
                        self.field_contribution_vector[site_index] @ dS_dtheta


        # Move in direction of force
        # If you want a physical interpretation, this is critically damped movement, where step_size_factor is dt
        # AKA 90s game engine physics, AKA Aristotlean physics, F = mv
        #
        # As we're in a coordinate system around current location, alpha = delta alpha

        alpha = step_size_factor * forces_free_alpha
        beta = step_size_factor * forces_free_beta
        theta = step_size_factor * forces_planar

        # print("alpha change:", alpha)
        # print("beta change:", beta)

        # Get the new moments
        ## Free sites

        new_moment_data = self.moment_data.copy()

        cos_beta = np.cos(beta)
        unrotated_moments = np.array([
            np.sin(beta),
            -np.sin(alpha) * cos_beta,
            np.cos(alpha) * cos_beta
        ])

        for component_index in range(self.n_components):
            for param_index, (site_index, _) in enumerate(self.free_sites):
                new_moment_data[site_index, component_index, :] = self.magnitudes[site_index, component_index] * \
                    rotation_matrices[param_index][component_index] @ unrotated_moments[:, param_index, component_index]

        ## Planar sites, apply the rotations to existing moments
        for component_index in range(self.n_components):
            for param_index, ((site_index, _), axis) in enumerate(zip(self.planar_sites, self.planar_axes)):
                new_moment_data[site_index, component_index, :] = \
                    rotation_matrix(theta[param_index, component_index], axis) \
                        @ self.moment_data[site_index, component_index, :]

        # Update
        self.moment_data = new_moment_data
