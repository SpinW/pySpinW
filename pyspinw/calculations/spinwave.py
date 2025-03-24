from dataclasses import dataclass

import numpy as np

from pyspinw.checks import check_sizes

@dataclass
class Coupling:
    index1: int
    index2: int
    matrix: np.ndarray


@check_sizes(rotations=('n_sites',3,3),
             magnitudes=('n_sites',),
             q_vectors=('n_q', 3),
             coupling_r=('n_sites', 3))
def spinwave_calculation(rotations: np.ndarray,
                         magnitudes: np.ndarray,
                         q_vectors: np.ndarray,
                         coupling_r: np.ndarray,
                         couplings: list[Coupling]):
    """ Main calculation step

    Unlike the main interface it takes indexed arrays, the meaning of the arrays is set elsewhere

    """

    n_sites = rotations.shape[0]

    # matrix used for calculating the phase factors, q.r in exp(i q.r) - could be too big, but we'll see
    qr = coupling_r.reshape(-1, 3) @ q_vectors.reshape(3, -1)

    phase_factors = np.exp(1j * qr)

    z = rotations[:,:,0] + 1j*rotations[:,:,1] # n-by-3, complex
    z_conj = z.conj()
    eta = rotations[:,:,2] # n-by-3 real

    # Create the A, B, and C matrices

    # Create the matrix sqrt(S_i S_j)/2
    root_mags = np.sqrt(0.5*magnitudes) # S_i / sqrt(2)
    total_spin_coefficients = root_mags.reshape(-1, 1) * root_mags.reshape(1, -1)

    A = np.zeros((n_sites, n_sites), dtype=complex)
    B = np.zeros((n_sites, n_sites), dtype=complex)
    C = np.zeros((n_sites, n_sites), dtype=complex)

    # Add the terms up to the total spin coefficient
    for coupling in couplings:
        i = coupling.index1
        j = coupling.index2
        A[i, j] += z[i, :] @ coupling.matrix @ z_conj[j, :]
        B[i, j] += z[i, :] @ coupling.matrix @ z[j, :]

        for l in range(n_sites):
            C[i, l] += eta[i, :] @ coupling.matrix @ eta[l, :]
            C[l, j] += eta[j, :] @ coupling.matrix @ eta[l, :]

    A *= total_spin_coefficients
    B *= total_spin_coefficients
    C *= total_spin_coefficients



if __name__ == "__main__":
    # Basic ferromagnet for checking
    # Single site
    rotations = np.eye(3).reshape(1,3,3)
    magnitudes = np.array([1.5]) # spin 3/2
    q_vectors = np.array([1,0,0]).reshape(1,3) * np.linspace(0, 10, 100).reshape(-1,1)
    coupling_r = np.array([0,0,0]).reshape(1,3)
    couplings = [Coupling(0,0,np.eye(3))]

    spinwave_calculation(rotations,
                         magnitudes,
                         q_vectors,
                         coupling_r,
                         couplings)


