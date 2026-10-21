""" Optimisation of propagation vectors """
import numpy as np

from pyspinw.hamiltonian import Hamiltonian
from pyspinw.symmetry.supercell import RotationSupercell


class PropagationVectorOptimisation:
    def __init__(self, hamiltonian: Hamiltonian):

        if not isinstance(hamiltonian.structure.supercell, RotationSupercell):
            raise TypeError("We can only optimise the propagation vector in incommensurate structures "
                            "as they are required to take continuous values")

        self.hamiltonian = hamiltonian
        self.supercell: RotationSupercell = hamiltonian.structure.supercell

        self.sites = hamiltonian.structure.sites
        self.n_sites = len(self.sites)

        self.site_uid_index_lookup = {site.unique_id: index for index, site in enumerate(self.sites)}

        #
        # Get basis vectors
        #
        basis_0 = self.supercell.propagation_vector.vector
        basis_0 /= np.sqrt(np.sum(basis_0**2))

        basis_1 = self.supercell.perpendicular
        basis_1 /= np.sqrt(np.sum(basis_1**2))

        basis_2 = np.cross(basis_0, basis_1)

        if not np.isclose(np.sum(basis_2**2), 1):
            raise ValueError("Perpendicular vector is not perpendicular to propagation vector")

        self.basis = np.array([basis_1, basis_2], dtype=complex)

    def minimise(self):
        pass


    def objective_function(self, vector_in_basis: np.ndarray):

        #
        # Build the coupling/anisotropy tensor
        #

        k_vector = self.basis @ vector_in_basis

        optimisation_tensor = np.zeros((self.n_sites, self.n_sites, 3, 3))
        for coupling in self.hamiltonian.couplings:
            index_1 = self.site_uid_index_lookup[coupling.site_1.unique_id]
            index_2 = self.site_uid_index_lookup[coupling.site_2.unique_id]

            k_projection = np.dot(coupling.lattice_vector(), k_vector)
            phase_factor = np.exp(2j*np.pi*k_projection)

            optimisation_tensor[index_1, index_2, :, :] += coupling.coupling_matrix * phase_factor
            optimisation_tensor[index_2, index_1, :, :] += coupling.coupling_matrix.T * phase_factor.conj()

        for anisotropy in self.hamiltonian.anisotropies:
            index = self.site_uid_index_lookup[anisotropy.site.unique_id]

            # phase_factor = 1.0    # Because distance is 0

            optimisation_tensor[index, index, :, :] += anisotropy.anisotropy_matrix

        # The next step is to tile into a 3n_sites-by-3n_sites matrix
        optimisation_matrix = optimisation_tensor.transpose(0, 2, 1, 3).reshape(3 * self.n_sites, 3 * self.n_sites)

        # The find the minimum real component of the eigenvectors
        return np.min(np.real(np.linalg.eigvals(optimisation_matrix)))