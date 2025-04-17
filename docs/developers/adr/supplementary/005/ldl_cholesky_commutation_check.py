""" Check that applying the commutation constraint works with an LDL decomposition """

from dataclasses import dataclass

import numpy as np

from numpy.linalg import cholesky
from scipy.linalg import ldl

from random_matrices import generate_random_hermitian_matrix

@dataclass
class Results:
    in_order: bool
    recovered_chol_matches: bool
    recovered_ldl_matches: bool
    decomposition_matches: bool
    diagonal_d: bool
    raw_eigenvalues_match: bool
    raw_eigenvectors_match: bool
    ordered_eigenvalues_match: bool
    ordered_eigenvectors_match: bool

    def failure_messages(self):
        messages = []
        if not self.recovered_chol_matches:
            messages.append("Recovered Cholesky not correct")
        if not self.recovered_ldl_matches:
            messages.append("Recovered LDL not correct")
        if not self.ordered_eigenvalues_match:
            messages.append("Ordered eigenvalues do not match")

        return messages


    def difference_messages(self):
        messages = []
        if not self.in_order:
            messages.append("Not in order")
        if not self.decomposition_matches:
            messages.append("Square root matrices not the same")
        if not self.raw_eigenvalues_match:
            messages.append("Raw eigenvalues do not match")
        if not self.raw_eigenvectors_match:
            messages.append("Raw eigenvectors do not match")
        if not self.ordered_eigenvectors_match:
            messages.append("Ordered eigenvectors do not match")

        return messages


    def status(self):
        failure_messages = self.failure_messages()
        if failure_messages:
            print("Non-match: " + ", ".join(self.difference_messages() + failure_messages))
            return False
        else:
            print("Match: " + ", ".join(self.difference_messages()))
            return True

def normalise(vectors):
    s = np.sqrt(np.sum(np.abs(vectors)**2, axis=1))
    return vectors / s.reshape(-1, 1)


def check_commutation_application(m):
    """ Check the commutation relationship """

    assert m.shape[0] == m.shape[1]
    assert m.shape[0] % 2 == 0

    n_states = m.shape[0] // 2

    delta = np.eye(2*n_states)
    delta[n_states:, :] *= -1

    # Cholesky version
    chol_sqrt_m = cholesky(m)
    recovered_chol = chol_sqrt_m @ chol_sqrt_m.conj().T
    recovered_chol_matches = np.all(np.abs(recovered_chol - m) < 1e-8)

    # Check 1: do we get the same decomposition with LDL
    l, d, perm = ldl(m)

    # Find out if the permutation is non-trivial
    diffs = perm[1:] - perm[:-1]
    in_order =  np.all(diffs == 1)

    # Find out if the 'diagonal' matrix contains block elements
    diagonal_d = np.all(np.abs(d - np.diag(np.diag(d))) < 1e-8)

    d_root = np.sqrt(d)

    ldl_sqrt_m = l @ d_root

    recovered_ldl = ldl_sqrt_m @ ldl_sqrt_m.conj().T
    recovered_ldl_matches = np.all(np.abs(recovered_ldl - m) < 1e-8)

    # Matrix to diagonalise
    kk_chol = chol_sqrt_m.conj().T @ delta @ chol_sqrt_m
    kk_ldl = ldl_sqrt_m.conj().T @ delta @ ldl_sqrt_m

    # Diagonalisation

    chol_eigenvalues, chol_eigenvectors = np.linalg.eig(kk_chol)
    ldl_eigenvalues, ldl_eigenvectors = np.linalg.eig(kk_ldl)

    # Check diagonalisation as is
    raw_eigenvalues_match = np.all(np.abs(chol_eigenvalues - ldl_eigenvalues) < 1e-8)
    raw_eigenvectors_match = np.all(np.abs(chol_eigenvectors - ldl_eigenvectors) < 1e-8)

    # Check ordered versions
    chol_ordering = np.argsort(chol_eigenvalues)
    ldl_ordering = np.argsort(ldl_eigenvalues)

    ordered_chol_eigenvalues = chol_eigenvalues[chol_ordering]
    ordered_ldl_eigenvalues = ldl_eigenvalues[ldl_ordering]
    ordered_eigenvalues_match = np.all(np.abs(ordered_chol_eigenvalues - ordered_ldl_eigenvalues) < 1e-8)

    ordered_chol_eigenvectors = normalise(chol_eigenvectors[chol_ordering, :])
    ordered_ldl_eigenvectors = normalise(ldl_eigenvectors[ldl_ordering, :])

    ordered_eigenvectors_match = np.all(np.abs(ordered_chol_eigenvectors - ordered_ldl_eigenvectors) < 1e-8)

    if not ordered_eigenvectors_match:
        print("ldl eigenvectors")
        print(ordered_ldl_eigenvectors)
        print("chol eigenvectors")
        print(ordered_chol_eigenvectors)

    # if False:
    #     print("Cholesky Real")
    #     print(chol_sqrt_m.real)
    #     print("LDL Real")
    #     print(ldl_sqrt_m.real)
    #     print("Cholesky Imag")
    #     print(chol_sqrt_m.imag)
    #     print("LDL Imag")
    #     print(ldl_sqrt_m.imag)
    #     print("Diff Real")
    #     print(np.abs(chol_sqrt_m.real - ldl_sqrt_m.real) > 1e-9)
    #     print("Diff Imag")
    #     print(np.abs(chol_sqrt_m.imag - ldl_sqrt_m.imag) > 1e-9)

    decomposition_matches = np.all(np.abs(chol_sqrt_m - ldl_sqrt_m) < 1e-8)

    return Results(
        recovered_chol_matches=recovered_chol_matches,
        recovered_ldl_matches=recovered_ldl_matches,
        in_order=in_order,
        decomposition_matches=decomposition_matches,
        diagonal_d=diagonal_d,
        raw_eigenvectors_match=raw_eigenvectors_match,
        raw_eigenvalues_match=raw_eigenvalues_match,
        ordered_eigenvalues_match=ordered_eigenvalues_match,
        ordered_eigenvectors_match=ordered_eigenvectors_match
    )

n = 1000
failed = 0
for i in range(n):
    m = generate_random_hermitian_matrix(2, seed=1000+i)

    results = check_commutation_application(m)

    success = results.status()

    if not success:
        failed += 1

print(f"{failed}/{n} did not match")

