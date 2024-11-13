""" Functions for working with different basis vectors"""

import numpy as np
from pyspinw.dimensionality import dimensionality_check

@dimensionality_check(vectors=(-1, 3))
def find_aligned_basis(vectors: np.ndarray, rcond: float | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """ Find a set of orthonormal basis vectors aligned with the first being aligned to the input vectors

    :param vectors: an n-by-3 numpy array containing the vectors to use to calculate the basis
    :param rcond: tolerance on the calculation, default is 3 * machine precision
    """

    if rcond is None:
        rcond = 3*np.finfo(vectors.dtype).eps # N * floating point epsilon

    # Lengths, for checking for zeros and for normalising
    lengths = np.sqrt(np.sum(vectors**2, axis=1))

    #
    # First basis vector will is just a normalised version of the input
    #

    # Avoid zero divisions by using np.divide instead of /
    zero_vectors = lengths < rcond
    e_1 = np.divide(vectors, lengths.reshape(-1, 1), where=~zero_vectors.reshape(-1, 1))

    # Assign zero vectors to z
    e_1[zero_vectors, 0] = 0.0
    e_1[zero_vectors, 1] = 0.0
    e_1[zero_vectors, 2] = 1.0

    #
    # Second basis, cross with x-axis vector, unless its pointing that way already, then we choose y-axis explicitly
    #

    x_aligned = np.abs(e_1[:, 0] - 1) < rcond

    x_vectors = np.zeros_like(vectors)
    x_vectors[:, 0] = 1.0

    e_2 = np.cross(e_1, x_vectors)
    e_2[x_aligned, 0] = 0.0
    e_2[x_aligned, 1] = 1.0
    e_2[x_aligned, 2] = 0.0
    e_2 /= np.sqrt(np.sum(e_2**2, axis=1)).reshape(-1, 1) # Normalise this one

    #
    # Third basis, just the cross of the first two
    #

    e_3 = np.cross(e_1, e_2) # Doesn't need normalising as e_1 and e_2 are orthogonal

    return e_1, e_2, e_3


def demo_find_aligned_basis():
    """ Example / test for find_aligned_basis"""

    test_vectors = np.array([
        [2, 0, 0],
        [0, 3, 0],
        [0, 0, 4],
        [5, 5, 5],
        [0, 0, 0]], dtype=float)

    e1, e2, e3 = find_aligned_basis(test_vectors)

    print(e1)
    print(e2)
    print(e3)


if __name__ == "__main__":
    demo_find_aligned_basis()
