"""Coupling Terms"""
from dataclasses import dataclass
from typing import ClassVar

import numpy as np

from pyspinw.cell_offsets import CellOffsetCoercible, CellOffset
from pyspinw.checks import check_sizes
from pyspinw.serialisation import SPWSerialisationContext, SPWSerialisable, numpy_serialise, SPWSerialisationError, \
    expects_keys, numpy_deserialise
from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.util import triple_product_matrix


@dataclass
class CouplingBaseDeserialisation:
    """ Class to hold basic properties of the coupling """
    name: str
    site_1: LatticeSite
    site_2: LatticeSite
    cell_offset: CellOffset


class Coupling(SPWSerialisable):
    """Coupling between different sites"""

    @check_sizes(coupling_matrix=(3,3))
    def __init__(self,
                 site_1: LatticeSite,
                 site_2: LatticeSite,
                 cell_offset: CellOffsetCoercible,
                 coupling_matrix: np.ndarray,
                 name: str = ""):

        self._name = name
        self._site_1 = site_1
        self._site_2 = site_2
        self._cell_offset = CellOffset.coerce(cell_offset)
        self._coupling_matrix = coupling_matrix


    coupling_type = "General"
    parameters: list[str] = []
    parameter_defaults: list[float] = []
    short_string = "M"


    @property
    def name(self):
        return self._name

    @property
    def site_1(self):
        return self._site_1

    @property
    def site_2(self):
        return self._site_2

    @property
    def cell_offset(self):
        return self._cell_offset

    @property
    def coupling_matrix(self) -> np.ndarray:
        """The coupling matrix for this coupling

        i.e. if H is the energy contribution for this coupling, S is the spin state, and
        M is the coupling matrix, we have

        H = S^T M S
        """
        return self._coupling_matrix

    @property
    def parameter_string(self) -> str:
        """ String representation of parameters """
        return ", ".join([f"{parameter}={self.__dict__[parameter]:.5g}" for parameter in self.parameters])

    @property
    def lattice_vector(self):
        """ Vector from site 1 to site 2 in lattice coordinates"""
        return self.cell_offset.vector + self._site_2.ijk - self._site_1.ijk

    def vector(self, unit_cell: UnitCell):
        """ Vector from site 1 to site 2 in cartesian coordinates (requires a unit cell definition)"""
        return unit_cell.fractional_to_cartesian(self.lattice_vector)

    def distance(self, unit_cell: UnitCell):
        """ Distance between sites """
        return np.sqrt(np.sum(self.vector(unit_cell)))

    def _base_serialisation(self, context: SPWSerialisationContext):
        return {
            "name": self.name,
            "site_1": self._site_1.serialise(context),
            "site_2": self._site_2.serialise(context),
            "cell_offset": self._cell_offset.serialise(context)
        }

    @staticmethod
    @expects_keys("name, site_1, site_2, cell_offset")
    def _base_deserialisation(json, context: SPWSerialisationContext):
        return CouplingBaseDeserialisation(
            name=json["name"],
            site_1=LatticeSite.deserialise(json["site_1"], context),
            site_2=LatticeSite.deserialise(json["site_2"], context),
            cell_offset=CellOffset.deserialise(json["cell_offset"], context))

    def serialise(self, context: SPWSerialisationContext) -> dict:
        return {
            "type": self.coupling_type.lower(),
            "data": self._coupling_serialise(context)
        }

    def _coupling_serialise(self, context: SPWSerialisationContext) -> dict:
        json = self._base_serialisation(context)
        json["matrix"] = numpy_serialise(self.coupling_matrix)
        return json

    @staticmethod
    @expects_keys("type, data")
    def deserialise(json: dict, context: SPWSerialisationContext):
        type_name = json["type"]
        return lowercase_coupling_lookup[type_name].deserialise(json["data"], context)


    @staticmethod
    @expects_keys("matrix")
    def _coupling_deserialise(data: dict, context: SPWSerialisationContext):
        base = Coupling._base_deserialisation(data, context)

        return Coupling(
            site_1=base.site_1,
            site_2=base.site_2,
            cell_offset=base.cell_offset,
            name=base.name,
            coupling_matrix=numpy_deserialise(data["matrix"]))



class HeisenbergCoupling(Coupling):
    """Heisenberg Coupling, which takes the form

    H_ij = J_ij (S_i . S_j)

    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param j: The coupling coefficient

    """

    coupling_type = "Heisenberg"
    parameters = ["j"]
    parameter_defaults = [1.0]
    short_string = "J"


    def __init__(self, site_1: LatticeSite, site_2: LatticeSite,
                 j: float,
                 cell_offset: CellOffsetCoercible = None,
                 name: str=""):

        self._j = j
        self._coupling_matrix = -j * np.eye(3)

        super().__init__(site_1=site_1,
                         site_2=site_2,
                         cell_offset=cell_offset,
                         coupling_matrix=self._coupling_matrix,
                         name=name)

    def _coupling_serialise(self, context: SPWSerialisationContext) -> dict:
        json = self._base_serialisation(context)
        json["j"] = self._j
        return json

    @staticmethod
    @expects_keys("j")
    def _coupling_deserialise(data: dict, context: SPWSerialisationContext):
        base = Coupling._base_deserialisation(data, context)

        return HeisenbergCoupling(
            site_1=base.site_1,
            site_2=base.site_2,
            cell_offset=base.cell_offset,
            name=base.name,
            j = data["j"])




class DiagonalCoupling(Coupling):
    """Diagonal coupling, which takes the form

    H_ij = Jxx_ij S^x_i S^x_j + Jyy_ij S^y_i S^y_j + Jzz_ij S^z_i S^z_j


    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param j: Vector containing the coupling coefficients for x, y, and z.

    """


    coupling_type = "Diagonal"
    parameters = ["j_x", "j_y", "j_z"]
    parameter_defaults = [1.0, 1.0, 1.0]
    short_string = "J"

    def __init__(self, site_1: LatticeSite, site_2: LatticeSite,
                 j_x: float, j_y: float, j_z: float,
                 cell_offset: CellOffsetCoercible = None,
                 name: str=""):

        self._j_x = j_x
        self._j_y = j_y
        self._j_z = j_z

        self._coupling_matrix = np.diag([j_x, j_y, j_z])

        super().__init__(site_1=site_1,
                         site_2=site_2,
                         cell_offset=cell_offset,
                         coupling_matrix=self._coupling_matrix,
                         name=name)

    def _coupling_serialise(self, context: SPWSerialisationContext) -> dict:
        json = self._base_serialisation(context)
        json["j_x"] = self._j_x
        json["j_y"] = self._j_y
        json["j_z"] = self._j_z
        return json

    @staticmethod
    @expects_keys("j_x, j_y, j_z")
    def _coupling_deserialise(data: dict, context: SPWSerialisationContext):
        base = Coupling._base_deserialisation(data, context)

        return DiagonalCoupling(
            site_1=base.site_1,
            site_2=base.site_2,
            cell_offset=base.cell_offset,
            name=base.name,
            j_x = data["j_x"],
            j_y = data["j_y"],
            j_z = data["j_z"])

class XYCoupling(Coupling):
    """ "XY"  coupling, which takes the form

    H_ij = J_ij (S^x_i S^x_j + S^y_i S^y_j)

    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param j: The coupling coefficient for the x and y components.

    """


    coupling_type = "XY"
    parameters = ["j"]
    parameter_defaults = [1.0]
    short_string = "J"

    def __init__(self, site_1: LatticeSite, site_2: LatticeSite,
                 j: float,
                 cell_offset: CellOffsetCoercible = None,
                 name: str=""):

        self._j = j
        self._coupling_matrix = np.eye(3)
        self._coupling_matrix[0,0] = -j
        self._coupling_matrix[1,1] = -j

        super().__init__(site_1=site_1,
                         site_2=site_2,
                         cell_offset=cell_offset,
                         coupling_matrix=self._coupling_matrix,
                         name=name)

    def _coupling_serialise(self, context: SPWSerialisationContext) -> dict:
        json = self._base_serialisation(context)
        json["j"] = self._j
        return json

    @staticmethod
    def _coupling_deserialise(data: dict, context: SPWSerialisationContext):
        base = Coupling._base_deserialisation(data, context)

        return XYCoupling(
            site_1=base.site_1,
            site_2=base.site_2,
            cell_offset=base.cell_offset,
            name=base.name,
            j = data["j"])



class XXZCoupling(Coupling):
    """ "XXZ" coupling, which takes the form

    H_ij = J_xy (S^x_i S^x_j + S^y_i S^y_j) + J_z (S^z_i S^z_j)

    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param j_xy: The coupling coefficient for the x and y components.
    :param j_z: The coupling coefficient for the z component

    """

    coupling_type = "XXZ"
    parameters = ["j_xy", "j_z"]
    parameter_defaults = [1.0, 1.0]
    short_string = "J"

    def __init__(self, site_1: LatticeSite, site_2: LatticeSite,
                 j_xy: float, j_z: float,
                 cell_offset: CellOffsetCoercible = None,
                 name: str=""):

        self._j_xy = j_xy
        self._j_z = j_z

        self._coupling_matrix = np.diag([j_xy, j_xy, j_z])

        super().__init__(site_1=site_1,
                         site_2=site_2,
                         cell_offset=cell_offset,
                         coupling_matrix=self._coupling_matrix,
                         name=name)

    def _coupling_serialise(self, context: SPWSerialisationContext) -> dict:
        json = self._base_serialisation(context)
        json["j_xy"] = self._j_xy
        json["j_z"] = self._j_z
        return json

    @staticmethod
    def _coupling_deserialise(data: dict, context: SPWSerialisationContext):
        base = Coupling._base_deserialisation(data, context)

        return XXZCoupling(
            site_1=base.site_1,
            site_2=base.site_2,
            cell_offset=base.cell_offset,
            name=base.name,
            j_xy = data["j_xy"],
            j_z = data["j_z"])

class IsingCoupling(Coupling):
    """Ising coupling (z component only), which takes the form

    H_ij = J_ij S^z_i S^z_j

    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param j: Scalar. The coupling coefficient J_ij.

    """

    coupling_type = "Ising"
    parameters = ["j_z"]
    parameter_defaults = [1.0]
    short_string = "J"

    def __init__(self, site_1: LatticeSite, site_2: LatticeSite,
                 j_z: float,
                 cell_offset: CellOffsetCoercible = None,
                 name: str=""):

        self._j_z = j_z

        self._coupling_matrix = np.diag([j_z])

        super().__init__(site_1=site_1,
                         site_2=site_2,
                         cell_offset=cell_offset,
                         coupling_matrix=self._coupling_matrix,
                         name=name)

    def _coupling_serialise(self, context: SPWSerialisationContext) -> dict:
        json = self._base_serialisation(context)
        json["j_z"] = self._j_z
        return json

    @staticmethod
    def _coupling_deserialise(data: dict, context: SPWSerialisationContext):
        base = Coupling._base_deserialisation(data, context)

        return IsingCoupling(
            site_1=base.site_1,
            site_2=base.site_2,
            cell_offset=base.cell_offset,
            name=base.name,
            j_z = data["j_z"])



class DMCoupling(Coupling):
    """Dzyaloshinskii–Moriya coupling, which takes the form

    H_ij = D_ij (S_i x S_j)

    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param d_x: x component of the d vector above
    :param d_y: x component of the d vector above
    :param d_z: x component of the d vector above

    """

    coupling_type = "Dzyaloshinskii–Moriya"
    parameters = ["d_x", "d_y", "d_z"]
    parameter_defaults = [1.0, 1.0, 1.0]
    short_string = "DM"


    def __init__(self, site_1: LatticeSite, site_2: LatticeSite,
                 d_x: float, d_y: float, d_z: float,
                 cell_offset: CellOffsetCoercible = None,
                 name: str=""):

        self._d_x = d_x
        self._d_y = d_y
        self._d_z = d_z

        self._coupling_matrix = triple_product_matrix(np.array([d_x, d_y, d_z], dtype=float))

        super().__init__(site_1=site_1,
                         site_2=site_2,
                         cell_offset=cell_offset,
                         coupling_matrix=self._coupling_matrix,
                         name=name)

    def _coupling_serialise(self, context: SPWSerialisationContext) -> dict:
        json = self._base_serialisation(context)
        json["d_x"] = self._d_x
        json["d_y"] = self._d_y
        json["d_z"] = self._d_z
        return json

    @staticmethod
    def _coupling_deserialise(data: dict, context: SPWSerialisationContext):
        base = Coupling._base_deserialisation(data, context)

        return DiagonalCoupling(
            site_1=base.site_1,
            site_2=base.site_2,
            cell_offset=base.cell_offset,
            name=base.name,
            j_x = data["d_x"],
            j_y = data["d_y"],
            j_z = data["d_z"])


couplings = [HeisenbergCoupling, DiagonalCoupling, XYCoupling, IsingCoupling, DMCoupling]
coupling_lookup = {coupling.coupling_type: coupling for coupling in couplings}
lowercase_coupling_lookup = {coupling.coupling_type.lower(): coupling for coupling in couplings}
