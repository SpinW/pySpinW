from dataclasses import dataclass

import numpy as np
from scipy.linalg import ldl

from pyspinw.checks import check_sizes

@dataclass
class Coupling:
    index1: int
    index2: int
    matrix: np.ndarray
    inter_site_vector: np.ndarray


@check_sizes(rotations=('n_sites',3,3),
             magnitudes=('n_sites',),
             q_vectors=('n_q', 3))
def spinwave_calculation(rotations: np.ndarray,
                         magnitudes: np.ndarray,
                         q_vectors: np.ndarray,
                         couplings: list[Coupling]):
    """ Main calculation step

    Unlike the main interface it takes indexed arrays, the meaning of the arrays is set elsewhere

    """

    n_sites = rotations.shape[0]

    z = rotations[:,:,0] + 1j*rotations[:,:,1] # n-by-3, complex
    z_conj = z.conj()
    eta = rotations[:,:,2] # n-by-3 real

    # Create the A, B, and C matrices

    # Create the matrix sqrt(S_i S_j)/2
    root_mags = np.sqrt(0.5*magnitudes) # S_i / sqrt(2)
    total_spin_coefficients = root_mags.reshape(-1, 1) * root_mags.reshape(1, -1)

    energies = []
    for q in q_vectors:

        A = np.zeros((n_sites, n_sites), dtype=complex)
        B = np.zeros((n_sites, n_sites), dtype=complex)
        C = np.zeros((n_sites, n_sites), dtype=complex)

        # Add the terms up to the total spin coefficient
        for coupling in couplings:

            phase_factor = np.exp((2j*np.pi)*np.dot(q, coupling.inter_site_vector))

            i = coupling.index1
            j = coupling.index2

            A[i, j] += z[i, :] @ coupling.matrix @ z_conj[j, :] * phase_factor
            B[i, j] += z[i, :] @ coupling.matrix @ z[j, :] * phase_factor


            for l in range(n_sites):
                C[i, l] += eta[i, :] @ coupling.matrix @ eta[l, :]
                C[l, j] += eta[j, :] @ coupling.matrix @ eta[l, :]

        A *= total_spin_coefficients
        B *= total_spin_coefficients
        C *= total_spin_coefficients # should be done over each site, but doesn't matter for now
        #
        # print("A", A)
        # print("B", B)
        # print("C", C)


        hamiltonian_matrix = np.block([[A - C, B], [B.conj().T, A - C]])

        #
        # We need to enforce the bosonic commutation properties, we do this
        # by finding the 'square root' of the matrix (i.e. finding K such that KK^dagger = H)
        # and then negating the second half.
        #
        # In matrix form we do
        #
        #     M = K C K^dagger
        #
        # where C is a diagonal matrix of length 2n, with the first n entries being 1, and the
        # remaining entries being -1.
        #
        # We can also do this via an LDL decomposition, but the method is slightly different
        # Instead we have a decomposition into

        try:
            sqrt_hamiltonian = np.linalg.cholesky(hamiltonian_matrix)

            sqrt_hamiltonian_with_commutation = sqrt_hamiltonian.copy()
            sqrt_hamiltonian_with_commutation[:, n_sites:] *= -1

            to_diagonalise = sqrt_hamiltonian @ np.conj(sqrt_hamiltonian_with_commutation).T

        except:
            l, d, perm = ldl(hamiltonian_matrix)

            to_negate = perm[n_sites:] # indexes


        eig_res = np.linalg.eig(to_diagonalise)

        energies.append(eig_res.eigenvalues) # These are currently the square energies

    return energies



if __name__ == "__main__":
    # Basic ferromagnet for checking
    # Single site
    rotations = np.eye(3).reshape(1,3,3)
    magnitudes = np.array([1.5]) # spin 3/2
    q_mags = np.linspace(0, 1, 100).reshape(-1,1)
    q_vectors = np.array([0,1,0]).reshape(1,3) * q_mags
    couplings = [Coupling(0,0,np.eye(3),inter_site_vector=np.array([0,1,0])),
                 Coupling(0, 0, np.eye(3), inter_site_vector=np.array([0, -1, 0])),
                 ]

    energies = spinwave_calculation(rotations,
                                     magnitudes,
                                     q_vectors,
                                     couplings)

    import matplotlib.pyplot as plt
    plt.plot(q_mags, energies)
    plt.show()
