""" Coupling Terms """

import numpy as np

from pyspinw.checks import check_sizes
from pyspinw._base import Coupling, Identifier
from pyspinw.util import triple_product_matrix


class HeisenbergCoupling(Coupling):
    """ Heisenberg Coupling, which takes the form

    H_ij = J_ij (S_i . S_j)

    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param j: The coupling coefficient

    """

    def __init__(self, site_1: Identifier, site_2: Identifier, j):
        super().__init__(site_1, site_2)

        self._j = j
        self._coupling_matrix = j * np.eye(3)


class DiagonalCoupling(Coupling):
    """ Diagonal coupling, which takes the form

    H_ij = J^x_ij S^x_i S^x_j + J^y_ij S^y_i S^y_j + J^z_ij S^z_i S^z_j


    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param j: Vector containing the coupling coefficients for x, y, and z.

    """

    @check_sizes(j=(3,))
    def __init__(self, site_1: Identifier, site_2: Identifier, j: np.ndarray):
        super().__init__(site_1, site_2)

        self._j = j
        self._coupling_matrix = np.diag(j)


class XYCoupling(Coupling):

    """ "XY" coupling, which takes the form

    H_ij = J_ij (S^x_i S^x_j + S^y_i S^y_j)

    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param j: The coupling coefficient for the x and y components.

    """

    def __init__(self, site_1: Identifier, site_2: Identifier, j):
        super().__init__(site_1, site_2)
        self._j = j
        self._coupling_matrix = np.diag([j, j, 0])



class XXZCoupling(Coupling):

    """ "XXZ" coupling, which takes the form

    H_ij = J_xy (S^x_i S^x_j + S^y_i S^y_j) + J_z (S^z_i S^z_j)

    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param j_xy: The coupling coefficient for the x and y components.
    :param j_z: The coupling coefficient for the z component

    """
    def __init__(self, site_1: Identifier, site_2: Identifier, j_xy, j_z):
        super().__init__(site_1, site_2)
        self._j_xy = j_xy
        self._j_z = j_z
        self._coupling_matrix = np.diag([j_xy, j_xy, j_z])


class IsingCoupling(Coupling):
    """ Ising coupling (z component only), which takes the form

    H_ij = J_ij S^z_i S^z_j

    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param j: Scalar. The coupling coefficient J_ij.

    """

    def __init__(self, site_1: Identifier, site_2: Identifier, j):
        super().__init__(site_1, site_2)
        self._j = j
        self._coupling_matrix = np.diag([0, 0, j])


class DMCoupling(Coupling):
    """ Dzyaloshinskii–Moriya coupling, which takes the form

    H_ij = D_ij (S_i x S_j)

    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param dm_vector: The vector D above

    """
    @check_sizes(d_vector=(3,), force_numpy=True)
    def __init__(self, site_1: Identifier, site_2: Identifier, dm_vector: np.ndarray):
        super().__init__(site_1, site_2)

        self._dm_vector = dm_vector
        self._coupling_matrix = triple_product_matrix(dm_vector)
