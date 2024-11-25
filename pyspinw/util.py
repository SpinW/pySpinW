""" Generally helpful functions that don't obviously live anywhere else in particular """

import numpy as np

from checks import check_sizes


@check_sizes(v=(3,), force_numpy=True)
def triple_product_matrix(v: np.ndarray):
    """Find the matrix, M, such that for all vectors X, Y

    X^T M Y = V . (X x Y)
    """

    x, y, z = v

    return np.array([
        [ 0,  z, -y],
        [-z,  0,  x],
        [ y, -x,  0]])


def demo_triple_product_matrix():
    """ Show an example of making a matrix that does the triple product """
    v = [1, 2, 3]
    m = triple_product_matrix(v)
    print(m)


if __name__ == "__main__":
    demo_triple_product_matrix()
