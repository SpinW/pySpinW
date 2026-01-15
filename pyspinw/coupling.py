"""Coupling Terms"""
from dataclasses import dataclass

import numpy as np

from pyspinw.cell_offsets import CellOffsetCoercible, CellOffset
from pyspinw.checks import check_sizes
from pyspinw.serialisation import SPWSerialisationContext, SPWSerialisable, numpy_serialise, \
    expects_keys, numpy_deserialise, SPWDeserialisationContext
from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.tolerances import tolerances
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

    serialisation_name = "coupling"

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
        """ Name of this coupling """
        return self._name

    @property
    def site_1(self):
        """ First site"""
        return self._site_1

    @property
    def site_2(self):
        """ Second site"""
        return self._site_2

    @property
    def cell_offset(self):
        """ Offset between unit cells containing sites 1 and 2"""
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
        # Note that we reference the _parameter value, not the property that references it
        substrings = []
        for parameter in self.parameters:
            value=self.__dict__["_" + parameter]
            substrings.append(f"{parameter}={value:.5g}")

        return ", ".join(substrings)

    def __repr__(self):

        direction_string = "<->" if self.is_symmetric() else "->"

        return "".join([
            self.__class__.__name__,
            f"('{self.name}', {self.site_1.name} {direction_string} {self.site_2.name}, offset={self.cell_offset}, ",
            self.parameter_string,
            ")"])

    @property
    def lattice_vector(self):
        """ Vector from site 1 to site 2 in lattice coordinates"""
        return self.cell_offset.vector + self._site_2.ijk - self._site_1.ijk

    def vector(self, unit_cell: UnitCell):
        """ Vector from site 1 to site 2 in cartesian coordinates (requires a unit cell definition)"""
        return unit_cell.fractional_to_cartesian(self.lattice_vector)

    def distance(self, unit_cell: UnitCell):
        """ Distance between sites """
        return np.sqrt(np.sum(self.vector(unit_cell)**2))

    def _base_serialisation(self, context: SPWSerialisationContext):
        return {
            "name": self.name,
            "site_1": self._site_1._serialise(context),
            "site_2": self._site_2._serialise(context),
            "cell_offset": self._cell_offset._serialise(context)
        }

    @check_sizes(coupling_matrix=(3,3), force_numpy=True, allow_nones=True)
    def updated(self,
                site_1 :LatticeSite | None = None,
                site_2: LatticeSite | None = None,
                cell_offset: CellOffset | None = None,
                name: str | None = None,
                coupling_matrix: np.ndarray | None = None):
        """ Get version of this site with specified parameters updated"""
        return Coupling(
            site_1=self.site_1 if site_1 is None else site_1,
            site_2=self.site_2 if site_2 is None else site_2,
            cell_offset=self.cell_offset if cell_offset is None else cell_offset,
            name=self.name if name is None else name,
            coupling_matrix=self.coupling_matrix if coupling_matrix is None else coupling_matrix)

    @staticmethod
    @expects_keys("name, site_1, site_2, cell_offset")
    def _base_deserialisation(json, context: SPWSerialisationContext):
        return CouplingBaseDeserialisation(
            name=json["name"],
            site_1=LatticeSite._deserialise(json["site_1"], context),
            site_2=LatticeSite._deserialise(json["site_2"], context),
            cell_offset=CellOffset._deserialise(json["cell_offset"], context))

    def _serialise(self, context: SPWSerialisationContext) -> dict:
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
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        type_name = json["type"]
        return lowercase_coupling_lookup[type_name]._coupling_deserialise(json["data"], context)


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

    def is_symmetric(self):
        """ Is this a symmetric coupling """
        return np.all(np.abs(self.coupling_matrix - self.coupling_matrix.T) < tolerances.IS_ZERO_TOL)


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
        self._coupling_matrix = j * np.eye(3)

        super().__init__(site_1=site_1,
                         site_2=site_2,
                         cell_offset=cell_offset,
                         coupling_matrix=self._coupling_matrix,
                         name=name)

    @property
    def j(self):
        """ Coupling constant """
        return self._j

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

    def is_symmetric(self):
        """ Is this a symmetric coupling """
        return True

    def updated(self,
                site_1 :LatticeSite | None = None,
                site_2: LatticeSite | None = None,
                cell_offset: CellOffset | None = None,
                name: str | None = None,
                j: float | None = None):
        """ Get version of this coupling with specified parameters updated"""
        return HeisenbergCoupling(
                site_1=self.site_1 if site_1 is None else site_1,
                site_2 = self.site_2 if site_2 is None else site_2,
                cell_offset = self.cell_offset if cell_offset is None else cell_offset,
                name = self.name if name is None else name,
                j = self.j if j is None else j)




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

    @property
    def j_x(self):
        """ Coupling constant for x """
        return self._j_x


    @property
    def j_y(self):
        """ Coupling constant for y """
        return self._j_y


    @property
    def j_z(self):
        """ Coupling constant for z """
        return self._j_z

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

    def is_symmetric(self):
        """ Is this a symmetric coupling """
        return True

    def updated(self,
                site_1 :LatticeSite | None = None,
                site_2: LatticeSite | None = None,
                cell_offset: CellOffset | None = None,
                name: str | None = None,
                j_x: float | None = None,
                j_y: float | None = None,
                j_z: float | None = None,
                ):
        """ Get version of this coupling with specified parameters updated"""
        return DiagonalCoupling(
                site_1=self.site_1 if site_1 is None else site_1,
                site_2 = self.site_2 if site_2 is None else site_2,
                cell_offset = self.cell_offset if cell_offset is None else cell_offset,
                name = self.name if name is None else name,
                j_x = self.j_x if j_x is None else j_x,
                j_y = self.j_y if j_y is None else j_y,
                j_z = self.j_z if j_z is None else j_z)

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


    @property
    def j(self):
        """ Coupling constant for x and y """
        return self._j

    def __init__(self, site_1: LatticeSite, site_2: LatticeSite,
                 j: float,
                 cell_offset: CellOffsetCoercible = None,
                 name: str=""):

        self._j = j
        self._coupling_matrix = np.diag([0.0, j, j], dtype=float)

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

    def is_symmetric(self):
        """ Is this a symmetric coupling """
        return True

    def updated(self,
                site_1: LatticeSite | None = None,
                site_2: LatticeSite | None = None,
                cell_offset: CellOffset | None = None,
                name: str | None = None,
                j: float | None = None):
        """ Get version of this coupling with specified parameters updated"""
        return XYCoupling(
                site_1=self.site_1 if site_1 is None else site_1,
                site_2=self.site_2 if site_2 is None else site_2,
                cell_offset=self.cell_offset if cell_offset is None else cell_offset,
                name=self.name if name is None else name,
                j=self.j if j is None else j)

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



    @property
    def j_xy(self):
        """ Coupling constant for x and y """
        return self._j_xy


    @property
    def j_z(self):
        """ Coupling constant for z """
        return self._j_z

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


    def updated(self,
                site_1: LatticeSite | None = None,
                site_2: LatticeSite | None = None,
                cell_offset: CellOffset | None = None,
                name: str | None = None,
                j_xy: float | None = None,
                j_z: float | None = None,
                ):
        """ Get version of this coupling with specified parameters updated"""
        return XXZCoupling(
                site_1=self.site_1 if site_1 is None else site_1,
                site_2=self.site_2 if site_2 is None else site_2,
                cell_offset=self.cell_offset if cell_offset is None else cell_offset,
                name=self.name if name is None else name,
                j_xy=self.j_xy if j_xy is None else j_xy,
                j_z=self.j_z if j_z is None else j_z)
    def is_symmetric(self):
        """ Is this a symmetric coupling """
        return True

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

    @property
    def j_z(self):
        """ Coupling constant for z """
        return self._j_z

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


    def updated(self,
                site_1: LatticeSite | None = None,
                site_2: LatticeSite | None = None,
                cell_offset: CellOffset | None = None,
                name: str | None = None,
                j_z: float | None = None):
        """ Get version of this coupling with specified parameters updated"""
        return IsingCoupling(
                site_1=self.site_1 if site_1 is None else site_1,
                site_2=self.site_2 if site_2 is None else site_2,
                cell_offset=self.cell_offset if cell_offset is None else cell_offset,
                name=self.name if name is None else name,
                j_z=self.j_z if j_z is None else j_z)

    def is_symmetric(self):
        """ Is this a symmetric coupling """
        return True

class DMCoupling(Coupling):
    """Dzyaloshinskii–Moriya coupling, which takes the form

    H_ij = D_ij (S_i x S_j)

    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param d_x: x component of the d vector above
    :param d_y: x component of the d vector above
    :param d_z: x component of the d vector above

    """

    coupling_type = "Dzyaloshinskii-Moriya"
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


    @property
    def d_x(self):
        """ DM coupling constant for x """
        return self._d_x

    @property
    def d_y(self):
        """ DM coupling constant for y """
        return self._d_y

    @property
    def d_z(self):
        """ DM coupling constant for z """
        return self._d_z

    def _coupling_serialise(self, context: SPWSerialisationContext) -> dict:
        json = self._base_serialisation(context)
        json["d_x"] = self._d_x
        json["d_y"] = self._d_y
        json["d_z"] = self._d_z
        return json

    @staticmethod
    def _coupling_deserialise(data: dict, context: SPWSerialisationContext):
        base = Coupling._base_deserialisation(data, context)

        return DMCoupling(
            site_1=base.site_1,
            site_2=base.site_2,
            cell_offset=base.cell_offset,
            name=base.name,
            d_x = data["d_x"],
            d_y = data["d_y"],
            d_z = data["d_z"])


    def updated(self,
                site_1: LatticeSite | None = None,
                site_2: LatticeSite | None = None,
                cell_offset: CellOffset | None = None,
                name: str | None = None,
                d_x: float | None = None,
                d_y: float | None = None,
                d_z: float | None = None):
        """ Get version of this coupling with specified parameters updated"""
        return DMCoupling(
                site_1=self.site_1 if site_1 is None else site_1,
                site_2=self.site_2 if site_2 is None else site_2,
                cell_offset=self.cell_offset if cell_offset is None else cell_offset,
                name=self.name if name is not None else name,
                d_x=self.d_x if d_x is None else d_x,
                d_y=self.d_y if d_y is None else d_y,
                d_z=self.d_z if d_z is None else d_z)


    def is_symmetric(self):
        """ Is this a symmetric coupling """
        return self.d_x == 0 and self.d_y == 0 and self.d_z == 0


couplings = [HeisenbergCoupling, DiagonalCoupling, XYCoupling, IsingCoupling, DMCoupling]
coupling_lookup = {coupling.coupling_type: coupling for coupling in couplings}
lowercase_coupling_lookup = {coupling.coupling_type.lower(): coupling for coupling in couplings}
