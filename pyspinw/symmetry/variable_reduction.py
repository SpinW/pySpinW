import numpy as np

def rref(A, tol=1e-12):
    """ Reduced row echelon form - not in scipy or numpy"""

    A = A.astype(float).copy()
    rows, cols = A.shape
    r = 0

    for c in range(cols):
        if r >= rows:
            break

        pivot = np.argmax(np.abs(A[r:, c])) + r

        if abs(A[pivot, c]) < tol:
            continue

        A[[r, pivot]] = A[[pivot, r]]

        A[r] /= A[r, c]

        for i in range(rows):
            if i != r:
                A[i] -= A[i, c] * A[r]

        r += 1

    A[np.abs(A) < tol] = 0
    return A

def exchange_matrix_constraints(transformations: list[np.ndarray]):
    """ Get the constraints on the exchange matrix based on symmetry transformations """

    parts = [np.kron(transformation, transformation) - np.eye(9) for transformation in transformations]

