""" Unit cells"""

import numpy as np

from ase.geometry.cell import cellpar_to_cell

from pyspinw.serialisation import SPWSerialisable, SPWSerialisationContext, SPWDeserialisationContext, expects_keys, \
    SPWSerialisationError, numpy_serialise, numpy_deserialise, vec3_serialise


class BadCellDefinition(Exception):
    """ Error for bad cell definition"""

class RawUnitCell(SPWSerialisable):
    """ Unit cell defined in terms of a matrix, its subclass `UnitCell` is constructed by lengths and angles"""

    _unit_cell_name = "raw"

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

    def _unit_cell_serialise(self):
        return {"xyz": numpy_serialise(self._xyz)}

    @staticmethod
    @expects_keys("xyz")
    def _unit_cell_deserialise(json):
        return RawUnitCell(numpy_deserialise(json["xyz"]))

    def _serialise(self, context: SPWSerialisationContext) -> dict:
        return {
            "type": self._unit_cell_name,
            "data": self._unit_cell_serialise()
        }

    @staticmethod
    @expects_keys("type, data")
    def _deserialise(json, context: SPWDeserialisationContext):

        try:
            return unit_cell_types[json["type"]]._unit_cell_deserialise(json["data"])

        except KeyError as ke:
            names = ", ".join([f"'{key}'" for key in unit_cell_types])
            raise SPWSerialisationError(f"Expected unit cell type to be one of {names}") from ke



class UnitCell(RawUnitCell):
    """ A Unit Cell Definition """

    _unit_cell_name = "abc"

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


    def _unit_cell_serialise(self):
        return {
            "a": self.a,
            "b": self.b,
            "c": self.c,
            "alpha": self.alpha,
            "beta": self.beta,
            "gamma": self.gamma,
            "ab_normal": vec3_serialise(*self.ab_normal),
            "direction": vec3_serialise(*self.direction)
        }

    @staticmethod
    @expects_keys("a, b, c, alpha, beta, gamma, ab_normal, direction")
    def _unit_cell_deserialise(json):
        return UnitCell(
            a=json["a"],
            b=json["b"],
            c=json["c"],
            alpha=json["alpha"],
            beta=json["beta"],
            gamma=json["gamma"],
            ab_normal=json["ab_normal"],
            direction=json["direction"]
        )


unit_cell_types = {cls._unit_cell_name: cls for cls in [RawUnitCell, UnitCell]}
