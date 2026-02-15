"""Spinwave Calculations"""

import multiprocessing
import traceback
from concurrent.futures import wait, ProcessPoolExecutor
from dataclasses import dataclass
from enum import Enum
from typing import Optional

import numpy as np
from scipy.linalg import ldl, solve

from pyspinw.checks import check_sizes
from pyspinw.constants import MU_B

# smallest energy not considered negligible (in meV)
ZERO_ENERGY_TOL = 1e-12

# Disable linting for bad variable names, because they should match the docs
# ruff: noqa: E741


@dataclass
class Coupling:
    """Temporary description of the coupling between atoms"""

    index1: int
    index2: int
    matrix: np.ndarray
    inter_site_vector: np.ndarray


@dataclass
class MagneticField:
    """Description of an external magnetic field."""

    vector: np.ndarray
    g_tensors: list[np.ndarray]


@dataclass
class SpinwaveResult:
    """Results from a spinwave calculation"""

    q_vectors: np.ndarray
    raw_energies: np.ndarray
    intensities: Optional[np.ndarray] = None
    sab: Optional[np.ndarray] = None

def _calc_q_independent(
        rotations: list[np.ndarray],
        magnitudes: np.ndarray,
        couplings: list[Coupling],
        field: MagneticField | None = None):
    """Calculate the q-independent matrices for the spinwave calculation.

    Returns the C matrix, the z-components of rotations, the spin coefficients matrix,
    and Zeeman term if relevant.
    """
    rotations = np.array(rotations)
    n_sites = rotations.shape[0]

    z = rotations[:, :, 0] + 1j * rotations[:, :, 1]  # n-by-3, complex
    eta = rotations[:, :, 2]  # n-by-3 real

    # Create the matrix sqrt(S_i S_j)/2
    root_mags = np.sqrt(0.5 * magnitudes)  # S_i / sqrt(2)
    spin_coefficients = np.outer(root_mags, root_mags)

    # calculate the C matrix for h(q), which is q-independent
    C = np.zeros((n_sites, n_sites), dtype=complex)
    for coupling in couplings:
        i, j = (coupling.index1, coupling.index2)
        C[j, j] += magnitudes[j] * eta[i, :].T @ coupling.matrix @ eta[j, :]

    # calculate the Zeeman term for the A matrix (A^z in Toth & Lake)
    if field is not None:
        Az = np.array([field.vector @ g_tensor @ eta[i, :] for i, g_tensor in enumerate(field.g_tensors)])
        Az *= -1 / 2 * MU_B
        Az = np.diag(Az)
    else:
        Az = None

    return C, z, spin_coefficients, Az


def _calc_sqrt_hamiltonian(
        q: np.ndarray,
        C: np.ndarray,
        n_sites: int,
        z: np.ndarray,
        spin_coefficients: np.ndarray,
        couplings: list[Coupling],
        Az: np.ndarray | None = None):
    """Calculate the square root of the 'Hamiltonian' h(q) for a single q-value."""
    A = np.zeros((n_sites, n_sites), dtype=complex)
    B = np.zeros((n_sites, n_sites), dtype=complex)
    z_conj = z.conj()

    # Add the terms up to the total spin coefficient
    for coupling in couplings:
        phase_factor = np.exp((2j * np.pi) * np.dot(q, coupling.inter_site_vector))

        i = coupling.index1
        j = coupling.index2

        A[i, j] += z[i, :] @ coupling.matrix @ z_conj[j, :] * phase_factor
        B[i, j] += z[i, :] @ coupling.matrix @ z[j, :] * phase_factor

    A *= spin_coefficients
    B *= spin_coefficients

    if Az is not None:
        A += Az

    hamiltonian_matrix = np.block([[A - C, B], [B.conj().T, A.conj().T - C]])
    hamiltonian_matrix += np.diag(np.array(range(hamiltonian_matrix.shape[0]))*1e-12)

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

    # assert np.sum(np.imag(np.linalg.eig(hamiltonian_matrix)[0])) < 1e-6, 'Non hermitian Hamiltonian!'

    try:
        sqrt_hamiltonian = np.linalg.cholesky(hamiltonian_matrix)
    except np.linalg.LinAlgError:  # Catch postive definiteness errors
        # l, d, perm = ldl(hamiltonian_matrix) # To LDL^\dagger (i.e. adjoint on right)
        # TODO: Check for actual diagonal (could potentially contain non-diagonal 2x2 blocks)
        l, d, p = ldl(hamiltonian_matrix)  # To LDL^\dagger (i.e. adjoint on right)
        sqrt_hamiltonian = l @ np.sqrt(d)

    return sqrt_hamiltonian


def _get_q_chunks(q_vectors: np.ndarray, n_proc: int):
    nq = int(np.floor(q_vectors.shape[0] / n_proc))
    return [q_vectors[i * nq : (i + 1) * nq] for i in range(n_proc - 1)] + [q_vectors[(n_proc - 1) * nq :]]


def energies(
        rotations: list[np.ndarray],
        magnitudes: np.ndarray,
        q_vectors: np.ndarray,
        couplings: list[Coupling],
        field: MagneticField | None = None) \
            -> np.ndarray:
    """Calculate the spinwave energies for a set of q-vectors.

    Unlike the main interface it takes indexed arrays, the meaning of the arrays is set elsewhere
    """
    C, z, spin_coefficients, Az = _calc_q_independent(rotations, magnitudes, couplings, field)
    n_sites = len(rotations)

    # Linear algebra routines in numpy are already parallelised and usually use 4 cores
    # for a single process, so we want to reduce contention by using fewer processes.
    n_proc = max(int(np.floor(multiprocessing.cpu_count() / 4)), 1)
    with ProcessPoolExecutor() as executor:
        q_calculations = [
            executor.submit(_calc_chunk_energies, q, C, n_sites, z, spin_coefficients, couplings, Az)
            for q in _get_q_chunks(q_vectors, n_proc)
        ]
    wait(q_calculations)
    energies = np.concat(tuple(future.result() for future in q_calculations))

    # return SpinwaveResult( q_vectors=q_vectors, raw_energies=energies, method=[])
    return energies


def _calc_chunk_energies(
        q_vectors: np.ndarray,
        C: np.ndarray,
        n_sites: int,
        z: np.ndarray,
        spin_coefficients: np.ndarray,
        couplings: list[Coupling],
        Az: np.ndarray | None = None):
    """Calculate the energies for a set of q-values."""
    energies = []
    for q in q_vectors:
        sqrt_hamiltonian = _calc_sqrt_hamiltonian(q, C, n_sites, z, spin_coefficients, couplings, Az)

        sqrt_hamiltonian_with_commutation = sqrt_hamiltonian.copy()
        sqrt_hamiltonian_with_commutation[n_sites:, :] *= -1  # This is C*K

        to_diagonalise = np.conj(sqrt_hamiltonian).T @ sqrt_hamiltonian_with_commutation

        energies.append(np.linalg.eigvals(to_diagonalise))

    return energies


def spinwave_calculation(
        rotations: list[np.ndarray],
        magnitudes: np.ndarray,
        q_vectors: np.ndarray,
        couplings: list[Coupling],
        positions: list[np.ndarray],
        rlu_to_cart: np.ndarray = np.eye(3),
        field: MagneticField | None = None,
        save_sab: bool = False) -> tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """Calculate the energies and spin-spin correlation for a set of q-vectors."""
    C, z, spin_coefficients, Az = _calc_q_independent(rotations, magnitudes, couplings, field)
    n_sites = len(rotations)

    # Linear algebra routines in numpy are already parallelised and usually use 4 cores
    # for a single process, so we want to reduce contention by using fewer processes.
    n_proc = max(int(np.floor(multiprocessing.cpu_count() / 4)), 1)
    with ProcessPoolExecutor() as executor:
        q_calculations = [
            executor.submit(
                _calc_chunk_spinwave, q, C, n_sites, z, spin_coefficients, couplings, positions, rlu_to_cart, Az, save_sab
            )
            for q in _get_q_chunks(q_vectors, n_proc)
        ]

    wait(q_calculations)

    results = []
    for future in q_calculations:
        try:
            results.append(future.result())
        except Exception as e:
            traceback.print_exception(type(e), e, e.__traceback__)

    energies = np.concat(tuple(result[0] for result in results))
    intensities = np.concat([result[1] for result in results])

    if save_sab:
        sab = [result[2] for result in results]
        # return SpinwaveResult(q_vectors, energies, intensities, sab)
        return energies, intensities, sab

    # return SpinwaveResult(q_vectors, energies, intensities)
    return energies, intensities

def _calc_chunk_spinwave(
        q_vectors: np.ndarray,
        C: np.ndarray,
        n_sites: int,
        z: np.ndarray,
        spin_coefficients: np.ndarray,
        couplings: list[Coupling],
        positions: list[np.ndarray],
        rlu_to_cart: np.ndarray,
        Az: np.ndarray | None = None,
        save_sab: bool = False) \
            -> tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """Calculate the energies and S'^alpha,beta for a chunk of q-values."""
    energies = []
    intensities = []
    sabs = []

    for q in q_vectors:
        sqrt_hamiltonian = _calc_sqrt_hamiltonian(q, C, n_sites, z, spin_coefficients, couplings, Az)

        sqrt_hamiltonian_with_commutation = sqrt_hamiltonian.copy()
        sqrt_hamiltonian_with_commutation[n_sites:, :] *= -1  # This is C*K

        to_diagonalise = np.conj(sqrt_hamiltonian).T @ sqrt_hamiltonian_with_commutation

        eigvals, eigvecs = np.linalg.eigh(to_diagonalise)
        energies.append(eigvals)

        ## calculate block matrices [ Y Z ; V W ] for S'^alpha,beta
        # first we get phase factor matrix where `phase_factors_matrix[i,j] = exp(i q (r_i - r_j))`
        phase_factors = np.array([np.exp(2j * np.pi * (q @ pos)) for pos in positions])
        phase_factors_matrix = np.outer(phase_factors, np.conj(phase_factors))

        coefficients = 2 * spin_coefficients * phase_factors_matrix

        # we store sab_blocks as a 3x3 array of arrays indexed over alpha, beta
        sab_blocks = np.zeros((3, 3), dtype=object)

        # note one can show:
        # Y*[alpha, beta] = Y[beta, alpha]
        # Z*[alpha, beta] = V[beta, alpha]
        # V*[alpha, beta] = Z[beta, alpha]
        # W*[alpha, beta] = W[beta, alpha]
        # thus we only need to calculate one triangle and can fill in the rest by conjugation
        for alpha in range(3):
            z_alphas = z[:, alpha]
            for beta in range(alpha + 1):
                z_betas = z[:, beta]

                # note V is conj(Z) and W is conj(Y) before we multiply by phase factors
                y_ab = np.outer(z_alphas, np.conj(z_betas))
                z_ab = np.outer(z_alphas, z_betas)
                v_ab = np.conj(z_ab) * coefficients
                w_ab = np.conj(y_ab) * coefficients

                y_ab *= coefficients
                z_ab *= coefficients

                sab_blocks[alpha, beta] = np.block([[y_ab, z_ab], [v_ab, w_ab]])

                if beta < alpha:
                    y_ba = np.conj(y_ab.T)
                    z_ba = np.conj(v_ab.T)
                    v_ba = np.conj(z_ab.T)
                    w_ba = np.conj(w_ab.T)

                    sab_blocks[beta, alpha] = np.block([[y_ba, z_ba], [v_ba, w_ba]])

        ## calculate transformation matrix for spin-spin correlation function
        # this is T = K^-1 U sqrt(E) where E is the diagonal 2 * n_sites matrix of eigenvalues
        # where the first n_sites entries are sqrt(eigval) and the remaining are sqrt(-eigval)
        # for the eigenvalues of the Hamiltonian
        sqrt_E = np.sqrt(np.abs(eigvals.copy()))
        sqrt_E[np.where(sqrt_E < ZERO_ENERGY_TOL)] = 0

        try:
            # rather than inverting K explicitly, calculate T by solving KT = U sqrt(E)
            T = solve(sqrt_hamiltonian.conj().T, eigvecs @ np.diag(sqrt_E))
            assert not np.isnan(T).any(), "singular matrix"
        except (AssertionError, np.linalg.LinAlgError):
            # if K is singular, then add a small amount to the diagonal.
            kk = sqrt_hamiltonian.conj().T
            for jj in range(kk.shape[0]):
                kk[jj, jj] += 1e-7
            if np.linalg.cond(kk) > 1e16:
                T = np.zeros((2 * n_sites, 2 * n_sites)) * np.nan
            else:
                T = solve(kk, eigvecs @ np.diag(sqrt_E))

        # Apply transformation matrix to S'^alpha,beta block matrices T*[VW;YZ]T
        # and then we just take the diagonal elements as that's all we need for
        # S'^alpha,beta(k, omega) at each eigenvalue
        # this is a 3x3x2N array indexed by [alpha, beta, omega]
        sab = np.array([[np.diag(T.conj().T @ sab_blocks[alpha, beta] @ T) for alpha in range(3)] for beta in range(3)])
        sab /= 2 * n_sites

        # take perpendicular component s_perp of sab
        if (q_mag := np.linalg.norm(q)) == 0:
            norm_q = np.array([0, 0, 0])
        else:
            norm_q = q @ rlu_to_cart
            norm_q /= np.sqrt(np.sum(norm_q**2))
        perp_factor = np.eye(3) - np.outer(norm_q, norm_q)
        # the einsum here performs elementwise multiplication and then sums over alpha, beta
        s_perp = np.einsum("ijk,ij->k", sab, perp_factor)
        intensities.append(s_perp)

    if save_sab:
        return energies, intensities, sabs

    return energies, intensities
