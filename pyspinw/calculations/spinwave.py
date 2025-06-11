"""Spinwave Calculations"""

import multiprocessing
from concurrent.futures import wait, ProcessPoolExecutor
from dataclasses import dataclass
from enum import Enum

import numpy as np
from scipy.linalg import ldl

from pyspinw.checks import check_sizes

PARALLEL = True


# Disable linting for bad variable names, because they should match the docs
# ruff: noqa: E741

@dataclass
class Coupling:
    """Temporary description of the coupling between atoms"""

    index1: int
    index2: int
    matrix: np.ndarray
    inter_site_vector: np.ndarray

class CalculationMethod(Enum):
    """Type of method used (for debugging purposes)"""

    CHOLESKY = 0
    LDL = 1

@dataclass
class SpinwaveResult:
    """Results from a spinwave calculation"""

    q_vectors: np.ndarray
    raw_energies: list[np.ndarray]
    method: list[CalculationMethod]


@check_sizes(magnitudes=('n_sites',),
             q_vectors=('n_q', 3))
def spinwave_calculation(rotations: list[np.ndarray],
                         magnitudes: np.ndarray,
                         q_vectors: np.ndarray,
                         couplings: list[Coupling]):
    """Main calculation step

    Unlike the main interface it takes indexed arrays, the meaning of the arrays is set elsewhere

    """
    rotations = np.array(rotations)
    n_sites = rotations.shape[0]

    z = rotations[:,:,0] + 1j*rotations[:,:,1] # n-by-3, complex
    eta = rotations[:,:,2] # n-by-3 real

    # Create the matrix sqrt(S_i S_j)/2
    root_mags = np.sqrt(0.5*magnitudes) # S_i / sqrt(2)
    spin_coefficients = root_mags.reshape(-1, 1) * root_mags.reshape(1, -1)

    # calculate the C matrix for h(q), which is q-independent
    C = np.zeros((n_sites, n_sites), dtype=complex)
    for coupling in couplings:
        i, j = (coupling.index1, coupling.index2)
        C[j, j] += spin_coefficients[j,j] * eta[i, :].T @ coupling.matrix @ eta[j, :]
    C *= 2

    if PARALLEL:
        # Linear algebra routines in numpy are already parallelised and usually use 4 cores
        # for a single process, so we want to reduce contention by using fewer processes.
        n_proc = int(np.floor(multiprocessing.cpu_count() / 4))
        with ProcessPoolExecutor() as executor:
            q_calculations = [executor.submit(_calc_q_chunks, q, C, n_sites, z, spin_coefficients, couplings)
                              for q in _get_q_chunks(q_vectors, n_proc)]
        wait(q_calculations)
        energies = np.concat(tuple(future.result() for future in q_calculations))
    else:
        energies = _calc_q_chunks(q_vectors, C, n_sites, z, spin_coefficients, couplings)

    #return SpinwaveResult( q_vectors=q_vectors, raw_energies=energies, method=[])
    return energies


def _get_q_chunks(q_vectors: np.ndarray, n_proc: int):
    nq = int(np.floor(q_vectors.shape[0] / n_proc))
    return [q_vectors[i*nq:(i+1)*nq] for i in range(n_proc - 1)] + [q_vectors[(n_proc - 1)*nq:]]


def _calc_q_chunks(q_vectors: np.ndarray,
                   C: np.ndarray,
                   n_sites: int,
                   z: np.ndarray,
                   spin_coefficients: np.ndarray,
                   couplings: list[Coupling]):
    """Calculate the energies for the 'Hamiltonian' h(q) for a set of q-values."""
    energies = []
    for q in q_vectors:
        A = np.zeros((n_sites, n_sites), dtype=complex)
        B = np.zeros((n_sites, n_sites), dtype=complex)
        z_conj = z.conj()
    
        # Add the terms up to the total spin coefficient
        for coupling in couplings:
    
            phase_factor = np.exp((2j*np.pi)*np.dot(q, coupling.inter_site_vector))
    
            i = coupling.index1
            j = coupling.index2
    
            A[i, j] += z[i, :] @ coupling.matrix @ z_conj[j, :] * phase_factor
            B[i, j] += z[i, :] @ coupling.matrix @ z[j, :] * phase_factor
    
        A *= spin_coefficients
        B *= spin_coefficients
        #
        # print("A", A)
        # print("B", B)
        # print("C", C)
    
    
        # hamiltonian_matrix = np.block([[A - C, B], [B.conj().T, A.conj().T - C]])
        hamiltonian_matrix = np.block([[A - C, B], [B.conj().T, A.conj().T - C]])
    
        #
        # We need to enforce the bosonic commutation properties, we do this
        # by finding the 'square root' of the matrix (i.e. finding K such that KK^dagger = H)
        # and then negating the second half.
        #
        # In matrix form we do
        #
        #     M = K^dagger C K
        #
        # where C is a diagonal matrix of length 2n, with the first n entries being 1, and the
        # remaining entries being -1. Note the adjoint is on the other side to the definition in
        # of the decomposition
        #
        # We can also do this via an LDL decomposition, but the method is very slightly different
    
    
        try:
            sqrt_hamiltonian = np.linalg.cholesky(hamiltonian_matrix)
            method = CalculationMethod.CHOLESKY
    
        except np.linalg.LinAlgError: # Catch postive definiteness errors
            # l, d, perm = ldl(hamiltonian_matrix) # To LDL^\dagger (i.e. adjoint on right)
            l, d, _ = ldl(hamiltonian_matrix) # To LDL^\dagger (i.e. adjoint on right)
            sqrt_hamiltonian = l @ np.sqrt(d)
    
            # TODO: Check for actual diagonal (could potentially contain non-diagonal 2x2 blocks)
            method = CalculationMethod.LDL
    
    
        sqrt_hamiltonian_with_commutation = sqrt_hamiltonian.copy()
        sqrt_hamiltonian_with_commutation[n_sites:, :] *= -1  # This is C*K
    
        to_diagonalise = np.conj(sqrt_hamiltonian).T @ sqrt_hamiltonian_with_commutation
    
        energies.append(np.linalg.eigvals(to_diagonalise))

    return energies
