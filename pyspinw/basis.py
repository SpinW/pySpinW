"""Functions for working with different basis vectors"""

import numpy as np
from pyspinw.checks import check_sizes
from pyspinw.tolerances import tolerances
from pyspinw.util import triple_product_matrix


@check_sizes(vectors=(-1, 3))
def find_aligned_basis(vectors: np.ndarray, rcond: float | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Find a set of orthonormal basis vectors aligned with the first being aligned to the input vectors

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
    e_1[zero_vectors, :] = np.array([[0.0,0.0,1.0]])

    #
    # Second basis, cross with x-axis vector, unless its pointing that way already, then we choose y-axis explicitly
    #

    x_aligned = np.abs(e_1[:, 0] - 1) < rcond

    x_vectors = np.zeros_like(vectors)
    x_vectors[:, 0] = 1.0

    e_2 = np.cross(e_1, x_vectors)
    e_2[x_aligned, :] = np.array([[0.0, 1.0, 0.0]])
    e_2 /= np.sqrt(np.sum(e_2**2, axis=1)).reshape(-1, 1) # Normalise this one

    #
    # Third basis, just the cross of the first two
    #

    e_3 = np.cross(e_1, e_2) # Doesn't need normalising as e_1 and e_2 are orthogonal

    return e_1, e_2, e_3

def site_rotations(vectors: np.ndarray, rcond: float | None = None):
    """ Return an aligned basis as an n-by-3-by-3 matrix """

    e1, e2, e3 = find_aligned_basis(vectors, rcond)

    return np.stack([e3, e2, e1], axis=2)

@check_sizes(axis=(3,), force_numpy=True)
def angle_axis_rotation_matrix(angle_rad, axis: np.ndarray):
    """ Get a matrix for rotation by `angle_rad` radians around the axis `axis`"""
    mag = np.sqrt(np.sum(axis**2))

    if mag < tolerances.IS_ZERO_TOL:
        raise ValueError("Rotation axis has very small or zero magnitude")

    axis = np.array(axis, dtype=float) / mag

    cos_theta = np.cos(angle_rad)
    sin_theta = np.sin(angle_rad)

    part1 = axis.reshape(-1, 1) * axis.reshape(1, -1) * (1 - cos_theta)
    part2 = cos_theta * np.eye(3)
    part3 = triple_product_matrix(axis * sin_theta)

    return part1 + part2 + part3

def demo_find_aligned_basis():
    """Example / test for find_aligned_basis"""
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

def demo_angle_axis_rotation_matrix():
    """ Show some plots demonstrating the angle-axis rotations"""
    import matplotlib.pyplot as plt

    v = np.array([1,0,0])

    xs = []
    ys = []

    for i in range(11):

        r = angle_axis_rotation_matrix(2*i*np.pi/11, [0,0,1])

        v2 = r @ v

        xs.append(v2[0])
        ys.append(v2[1])

    plt.scatter(xs, ys)

    xs = []
    ys = []

    for i in range(11):

        r = angle_axis_rotation_matrix(2*i*np.pi/11, [1,1,1])

        v2 = r @ v

        xs.append(v2[0])
        ys.append(v2[1])


    plt.scatter(xs, ys)
    plt.show()




if __name__ == "__main__":
    demo_find_aligned_basis()
    demo_angle_axis_rotation_matrix()
