"""Coupling Terms"""
from typing import ClassVar

import numpy as np

from pyspinw.checks import check_sizes
from pyspinw._base import Coupling, Identifier
from pyspinw.util import triple_product_matrix


class HeisenbergCoupling(Coupling):
    """Heisenberg Coupling, which takes the form

    H_ij = J_ij (S_i . S_j)

    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param j: The coupling coefficient

    """

    j: float

    coupling_type: ClassVar[str] = "Heisenberg"
    parameters: ClassVar[list[str]] = ["j"]

    def model_post_init(self, __context):
        self._coupling_matrix = self.j * np.eye(3)


class DiagonalCoupling(Coupling):
    """Diagonal coupling, which takes the form

    H_ij = Jxx_ij S^x_i S^x_j + Jyy_ij S^y_i S^y_j + Jzz_ij S^z_i S^z_j


    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param j: Vector containing the coupling coefficients for x, y, and z.

    """

    j_x: float
    j_y: float
    j_z: float


    coupling_type: ClassVar[str] = "Diagonal"
    parameters: ClassVar[list[str]] = ["j_x", "j_y", "j_z"]

    def model_post_init(self, __context):
        self._coupling_matrix = np.diag([self.j_x, self.j_y, self.j_z])


class XYCoupling(Coupling):
    """ "XY"  coupling, which takes the form

    H_ij = J_ij (S^x_i S^x_j + S^y_i S^y_j)

    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param j: The coupling coefficient for the x and y components.

    """

    j: float


    coupling_type: ClassVar[str] = "XY"
    parameters: ClassVar[list[str]] = ["j"]

    def model_post_init(self, __context):
        self._coupling_matrix = np.diag([self.j, self.j, 0.0], dtype=float)



class XXZCoupling(Coupling):
    """ "XXZ" coupling, which takes the form

    H_ij = J_xy (S^x_i S^x_j + S^y_i S^y_j) + J_z (S^z_i S^z_j)

    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param j_xy: The coupling coefficient for the x and y components.
    :param j_z: The coupling coefficient for the z component

    """

    j_xy: float
    j_z: float

    coupling_type: ClassVar[str] = "XXZ"
    parameters: ClassVar[list[str]] = ["j_xy", "j_z"]

    def model_post_init(self, __context):
        self._coupling_matrix = np.diag([self.j_xy, self.j_xy, self.j_z])


class IsingCoupling(Coupling):
    """Ising coupling (z component only), which takes the form

    H_ij = J_ij S^z_i S^z_j

    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param j: Scalar. The coupling coefficient J_ij.

    """

    j_z: float

    coupling_type: ClassVar[str] = "Ising"

    parameters: ClassVar[list[str]] = ["j_z"]

    def model_post_init(self, __context):
        self._coupling_matrix = np.diag([0, 0, self.j_z])




class DMCoupling(Coupling):
    """Dzyaloshinskii–Moriya coupling, which takes the form

    H_ij = D_ij (S_i x S_j)

    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param d_x: x component of the d vector above
    :param d_y: x component of the d vector above
    :param d_z: x component of the d vector above

    """

    d_x: float
    d_y: float
    d_z: float

    coupling_type: ClassVar[str] = "Dzyaloshinskii–Moriya"
    parameters: ClassVar[list[str]] = ["d_x", "d_y", "d_z"]

    def model_post_init(self, __context):

        self._coupling_matrix = triple_product_matrix(np.array([self.d_x, self.d_y, self.d_z]))

    @property
    def parameter_string(self):
        parts = []
        for parameter in self.parameters:
            parts.append(parameter + repr(self.__dict__[parameter]))

        return ", ".join(parts)


couplings = [HeisenbergCoupling, DiagonalCoupling, XYCoupling, IsingCoupling, DMCoupling]