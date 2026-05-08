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
        self._xyz = xyz # Transformation from lattice units to cartesian

        # Inverse transformation

        try:
            self._xyz_inv = np.linalg.inv(xyz)

        except np.linalg.LinAlgError as e:
            raise BadCellDefinition(f"{self._xyz} is not invertible")

        #
        # Matrices to convert spins
        #

        a_vector = xyz[0, :]
        b_vector = xyz[1, :]
        c_vector = xyz[2, :]

        # The fist coordinate only needs to be normalised
        spin_a = a_vector / np.sqrt(np.sum(a_vector**2))

        # The second component is "in the ab plane and perpendicular to b"
        #  we want to take away the component of b that is in the direction of a

        b_in_a = np.dot(spin_a, b_vector)    # scalar projection of b into a

        spin_b = b_vector - b_in_a * spin_a  # unnormalised
        spin_b /= np.sqrt(np.sum(spin_b**2)) # normalise

        # assert np.dot(spin_a, spin_b) < 1e-10, "spin_a and spin_b should be orthogonal"

        # last one should just be the cross product, though we need to be careful about handedness
        spin_c = np.cross(spin_a, spin_b)

        # assert np.dot(c_vector, spin_c) > 0, "spin c vector should in the direction of input c"

        # build the spin matrix, the row/column format has been checked and is the same as xyz
        self._xyz_spins = np.array([spin_a, spin_b, spin_c])

        try:
            self._xyz_spins_inv = np.linalg.inv(self._xyz_spins)

        except np.linalg.LinAlgError as e:

            raise BadCellDefinition(f"{self._xyz} doesn't allow an invertible spin definition")

    # @check_sizes(points=(-1, 3))
    def lattice_units_to_cartesian(self, points: np.ndarray):
        """ Convert a list of points  from the fractional (ijk) type to cartesian (xyz) """
        return points @ self._xyz

    # @check_sizes(points=(-1, 3))
    def cartesian_to_lattice_units(self, points: np.ndarray):
        """ Convert a list of points from cartesian (xyz) to fractional (ijk) """
        return points @ self._xyz_inv

    def spin_lattice_units_to_cartesian(self, spins: np.ndarray):
        """ Convert a list of spins from lattice units to cartesian """
        return spins @ self._xyz_spins

    def spin_cartesian_to_lattice_units(self, spins: np.ndarray):
        """ Convert a list of spins from cartesian to lattice units """
        return spins @ self._xyz_spins_inv

    @property
    def centre(self):
        """ Point at centre of cell in cartesian coordinates"""
        return self.lattice_units_to_cartesian(np.array([[0.5, 0.5, 0.5]]))[0, :]

    @property
    def main_diagonal_length(self):
        """ Length of primary diagonal """
        return np.sqrt(np.sum(self.lattice_units_to_cartesian(np.array([[1, 1, 1]]))[0, :] ** 2))

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

    @property
    def text_summary(self):
        """ Formatted text """
        return f"({self.a}, {self.b}, {self.c}), alpha={self.alpha}, beta={self.beta}, gamma={self.gamma}"

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

if __name__ == "__main__":
    UnitCell(a=1, b=2, c=3, beta=40)
