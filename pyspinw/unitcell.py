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

    @check_sizes(points=(-1, 3))
    def fractional_to_cartesian(self, points: np.ndarray):
        """ Convert a list of points  from the fractional (ijk) type to cartesian (xyz) """

        return points @ self._xyz

    @check_sizes(points=(-1, 3))
    def cartesian_to_fractional(self, points: np.ndarray):
        """ Convert a list of points from cartesian (xyz) to  fractional (ijk) """

        return points @ self._xyz_inv

    @property
    def centre(self):
        """ Point at centre of cell in cartesian coordinates"""
        return self.fractional_to_cartesian(np.array([[0.5,0.5,0.5]]))[0, :]


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

        """

        See `ase.geometry.cell.cellpar_to_cell` for details of parameters

        """

        self.a = a
        self.b = b
        self.c = c

        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma

        self.abc = np.array([a,b,c])

        xyz = cellpar_to_cell([a,b,c,alpha,beta,gamma], ab_normal=ab_normal, a_direction=direction)

        super().__init__(xyz)

    def __repr__(self):
        return f"UnitCell({self.a}, {self.b}, {self.c} | {self.alpha}, {self.beta}, {self.gamma})"