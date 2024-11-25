from idlelib.squeezer import Squeezer

import numpy as np
from numpy.random.mtrand import Sequence

from checks import check_sizes
from pyspinw._base import Coupling, Identifier

class ScalarCoupling(Coupling):

    def __init__(self, site_1: Identifier, site_2: Identifier, j):
        super().__init__(site_1, site_2)

        self._j = j
        self._coupling_matrix = j * np.eye(3)

class DiagonalCoupling(Coupling):

    @check_sizes(j=(3,))
    def __init__(self, site_1: Identifier, site_2: Identifier, j: np.ndarray):
        super().__init__(site_1, site_2)

        self._j = j
        self._coupling_matrix = np.diag(j)


class XYCoupling(Coupling):
    def __init__(self, site_1: Identifier, site_2: Identifier, j):
        super().__init__(site_1, site_2)
        self._j = j
        self._coupling_matrix = np.diag([j, j, 0])



class XXZCoupling(Coupling):
    def __init__(self, site_1: Identifier, site_2: Identifier, j_xy, j_z):
        super().__init__(site_1, site_2)
        self._j_xy = j_xy
        self._j_z = j_z
        self._coupling_matrix = np.diag([j_xy, j_xy, j_z])


class IsingCoupling(Coupling):
    def __init__(self, site_1: Identifier, site_2: Identifier, j):
        super().__init__(site_1, site_2)
        self._j = j
        self._coupling_matrix = np.diag([0, 0, j])


class DMCoupling(Coupling):
    """ Coupling between two"""
    @check_sizes(d_vector=(3,), force_numpy=True)
    def __init__(self, site_1: Identifier, site_2: Identifier, dm_vector: np.ndarray):
        super().__init__(site_1, site_2)

        self._dm_vector = dm_vector