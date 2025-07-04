""" Unit cells"""

import numpy as np

from pyspinw.checks import check_sizes

from ase.geometry.cell import cellpar_to_cell

class BadCellDefinition(Exception):
    """ Error for bad cell defintion"""

class RawUnitCell:
    """ Unit cell defined in terms of a matrix, its subclass `UnitCell` is constructed by lengths and angles"""

    def __init__(self, xyz):
        self._xyz = xyz

        try:
            self._xyz_inv = np.linalg.inv(xyz)

        except np.linalg.LinAlgError as e:
            raise BadCellDefinition(f"{self._xyz}")

    # @check_sizes(points=(-1, 3))
    def fractional_to_cartesian(self, points: np.ndarray):
        """ Convert a list of points  from the fractional (ijk) type to cartesian (xyz) """
        return points @ self._xyz

    # @check_sizes(points=(-1, 3))
    def cartesian_to_fractional(self, points: np.ndarray):
        """ Convert a list of points from cartesian (xyz) to  fractional (ijk) """
        return points @ self._xyz_inv

    @property
    def centre(self):
        """ Point at centre of cell in cartesian coordinates"""
        return self.fractional_to_cartesian(np.array([[0.5,0.5,0.5]]))[0, :]

    @property
    def main_diagonal_length(self):
        """ Length of primary diagonal """
        return np.sqrt(np.sum(self.fractional_to_cartesian(np.array([[1,1,1]]))[0, :]**2))

    def __eq__(self, other: "UnitCell"):
        """Equality (approximate)"""
        return np.all(np.abs(self._xyz - other._xyz) < 1e-10)


class UnitCell(RawUnitCell):
    """ A Unit Cell Definition """

    def __init__(self,
                 a: float,
                 b: float,
                 c: float,
                 alpha: float = 90,
                 beta: float = 90,
                 gamma: float = 90,
                 ab_normal: tuple[float, float, float]=(0,0,1),
                 direction: tuple[float, float, float] | None=None):
        """See `ase.geometry.cell.cellpar_to_cell` for details of parameters"""
        self.a = a
        self.b = b
        self.c = c

        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

        self.ab_normal = ab_normal
        self.direction = direction

        self.abc = np.array([a,b,c])

        xyz = cellpar_to_cell([a,b,c,alpha,beta,gamma], ab_normal=ab_normal, a_direction=direction)

        super().__init__(xyz)

    def updated(self,
                a: float | None = None,
                b: float | None = None,
                c: float | None = None,
                alpha: float | None = None,
                beta: float | None = None,
                gamma: float | None = None,
                ab_normal: tuple[float, float, float] | None = None,
                direction: tuple[float, float, float] | None = None,
                replace_direction=False):
        """ Create a new unit cell with updated parameters

        :returns: a new unit cell

        """
        a = self.a if a is None else a
        b = self.b if b is None else b
        c = self.c if c is None else c

        alpha = self.alpha if alpha is None else alpha
        beta = self.beta if beta is None else beta
        gamma = self.gamma if gamma is None else gamma

        ab_normal = self.ab_normal if ab_normal is None else ab_normal

        if not replace_direction:
            direction = self.direction

        return UnitCell(a, b, c, alpha, beta, gamma, ab_normal, direction)


    def __repr__(self):
        """repr implementation"""
        return f"UnitCell({self.a}, {self.b}, {self.c} | {self.alpha}, {self.beta}, {self.gamma})"
