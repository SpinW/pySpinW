"""Exchange Terms"""
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

_exchange_id_counter = -1
def _generate_unique_exchange_id():
    """ Generate a unique ID for each site currently loaded"""
    global _exchange_id_counter # noqa: PLW0603
    _exchange_id_counter += 1
    return _exchange_id_counter

@dataclass
class ExchangeBaseDeserialisation:
    """ Class to hold basic properties of the exchange """

    name: str
    site_1: LatticeSite
    site_2: LatticeSite
    cell_offset: CellOffset


class Exchange(SPWSerialisable):
    """Exchange between different sites"""

    serialisation_name = "exchange"

    @check_sizes(exchange_matrix=(3,3))
    def __init__(self,
                 site_1: LatticeSite,
                 site_2: LatticeSite,
                 cell_offset: CellOffsetCoercible,
                 exchange_matrix: np.ndarray,
                 name: str = ""):

        self._name = name
        self._site_1 = site_1
        self._site_2 = site_2
        self._cell_offset = CellOffset.coerce(cell_offset)
        self._exchange_matrix = exchange_matrix

        self._unique_id = _generate_unique_exchange_id()


    exchange_type = "General"
    parameters: list[str] = []
    parameter_defaults: list[float] = []
    short_string = "M"


    @property
    def name(self):
        """ Name of this exchange """
        return self._name

    @property
    def unique_id(self) -> int:
        """ Unique identifier for this exchange """
        return self._unique_id

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
    def exchange_matrix(self) -> np.ndarray:
        """The exchange matrix for this exchange

        i.e. if H is the energy contribution for this exchange, S is the spin state, and
        M is the exchange matrix, we have

        H = S^T M S
        """
        return self._exchange_matrix

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

    @check_sizes(exchange_matrix=(3,3), force_numpy=True, allow_nones=True)
    def updated(self,
                site_1 :LatticeSite | None = None,
                site_2: LatticeSite | None = None,
                cell_offset: CellOffset | None = None,
                name: str | None = None,
                exchange_matrix: np.ndarray | None = None):
        """ Get version of this site with specified parameters updated"""
        return Exchange(
            site_1=self.site_1 if site_1 is None else site_1,
            site_2=self.site_2 if site_2 is None else site_2,
            cell_offset=self.cell_offset if cell_offset is None else cell_offset,
            name=self.name if name is None else name,
            exchange_matrix=self.exchange_matrix if exchange_matrix is None else exchange_matrix)

    @staticmethod
    @expects_keys("name, site_1, site_2, cell_offset")
    def _base_deserialisation(json, context: SPWSerialisationContext):
        return ExchangeBaseDeserialisation(
            name=json["name"],
            site_1=LatticeSite._deserialise(json["site_1"], context),
            site_2=LatticeSite._deserialise(json["site_2"], context),
            cell_offset=CellOffset._deserialise(json["cell_offset"], context))

    def _serialise(self, context: SPWSerialisationContext) -> dict:
        return {
            "type": self.exchange_type.lower(),
            "data": self._exchange_serialise(context)
        }

    def _exchange_serialise(self, context: SPWSerialisationContext) -> dict:
        json = self._base_serialisation(context)
        json["matrix"] = numpy_serialise(self._exchange_matrix)
        return json

    @staticmethod
    @expects_keys("type, data")
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        type_name = json["type"]
        return lowercase_exchange_lookup[type_name]._exchange_deserialise(json["data"], context)


    @staticmethod
    @expects_keys("matrix")
    def _exchange_deserialise(data: dict, context: SPWSerialisationContext):
        base = Exchange._base_deserialisation(data, context)

        return Exchange(
            site_1=base.site_1,
            site_2=base.site_2,
            cell_offset=base.cell_offset,
            name=base.name,
            exchange_matrix=numpy_deserialise(data["matrix"]))

    def is_symmetric(self):
        """ Is this a symmetric exchange """
        return np.all(np.abs(self.exchange_matrix - self.exchange_matrix.T) < tolerances.IS_ZERO_TOL)


class HeisenbergExchange(Exchange):
    """Heisenberg Exchange, which takes the form

    H_ij = J_ij (S_i . S_j)

    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param j: The exchange coefficient

    """

    exchange_type = "Heisenberg"
    parameters = ["j"]
    parameter_defaults = [1.0]
    short_string = "J"


    def __init__(self, site_1: LatticeSite, site_2: LatticeSite,
                 j: float,
                 cell_offset: CellOffsetCoercible = None,
                 name: str=""):

        self._j = j
        self._exchange_matrix = j * np.eye(3)

        super().__init__(site_1=site_1,
                         site_2=site_2,
                         cell_offset=cell_offset,
                         exchange_matrix=self._exchange_matrix,
                         name=name)

    @property
    def j(self):
        """ Exchange constant """
        return self._j

    def _exchange_serialise(self, context: SPWSerialisationContext) -> dict:
        json = self._base_serialisation(context)
        json["j"] = self._j
        return json

    @staticmethod
    @expects_keys("j")
    def _exchange_deserialise(data: dict, context: SPWSerialisationContext):
        base = Exchange._base_deserialisation(data, context)

        return HeisenbergExchange(
            site_1=base.site_1,
            site_2=base.site_2,
            cell_offset=base.cell_offset,
            name=base.name,
            j = data["j"])

    def is_symmetric(self):
        """ Is this a symmetric exchange """
        return True

    def updated(self,
                site_1 :LatticeSite | None = None,
                site_2: LatticeSite | None = None,
                cell_offset: CellOffset | None = None,
                name: str | None = None,
                j: float | None = None):
        """ Get version of this exchange with specified parameters updated"""
        return HeisenbergExchange(
                site_1=self.site_1 if site_1 is None else site_1,
                site_2 = self.site_2 if site_2 is None else site_2,
                cell_offset = self.cell_offset if cell_offset is None else cell_offset,
                name = self.name if name is None else name,
                j = self.j if j is None else j)




class DiagonalExchange(Exchange):
    """Diagonal exchange, which takes the form

    H_ij = Jxx_ij S^x_i S^x_j + Jyy_ij S^y_i S^y_j + Jzz_ij S^z_i S^z_j


    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param j: Vector containing the exchange coefficients for x, y, and z.

    """

    exchange_type = "Diagonal"
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

        self._exchange_matrix = np.diag([j_x, j_y, j_z])

        super().__init__(site_1=site_1,
                         site_2=site_2,
                         cell_offset=cell_offset,
                         exchange_matrix=self._exchange_matrix,
                         name=name)

    @property
    def j_x(self):
        """ Exchange constant for x """
        return self._j_x


    @property
    def j_y(self):
        """ Exchange constant for y """
        return self._j_y


    @property
    def j_z(self):
        """ Exchange constant for z """
        return self._j_z

    def _exchange_serialise(self, context: SPWSerialisationContext) -> dict:
        json = self._base_serialisation(context)
        json["j_x"] = self._j_x
        json["j_y"] = self._j_y
        json["j_z"] = self._j_z
        return json

    @staticmethod
    @expects_keys("j_x, j_y, j_z")
    def _exchange_deserialise(data: dict, context: SPWSerialisationContext):
        base = Exchange._base_deserialisation(data, context)

        return DiagonalExchange(
            site_1=base.site_1,
            site_2=base.site_2,
            cell_offset=base.cell_offset,
            name=base.name,
            j_x = data["j_x"],
            j_y = data["j_y"],
            j_z = data["j_z"])

    def is_symmetric(self):
        """ Is this a symmetric exchange """
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
        """ Get version of this exchange with specified parameters updated"""
        return DiagonalExchange(
                site_1=self.site_1 if site_1 is None else site_1,
                site_2 = self.site_2 if site_2 is None else site_2,
                cell_offset = self.cell_offset if cell_offset is None else cell_offset,
                name = self.name if name is None else name,
                j_x = self.j_x if j_x is None else j_x,
                j_y = self.j_y if j_y is None else j_y,
                j_z = self.j_z if j_z is None else j_z)

class XYExchange(Exchange):
    """ "XY"  exchange, which takes the form

    H_ij = J_ij (S^x_i S^x_j + S^y_i S^y_j)

    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param j: The exchange coefficient for the x and y components.

    """

    exchange_type = "XY"
    parameters = ["j"]
    parameter_defaults = [1.0]
    short_string = "J"


    @property
    def j(self):
        """ Exchange constant for x and y """
        return self._j

    def __init__(self, site_1: LatticeSite, site_2: LatticeSite,
                 j: float,
                 cell_offset: CellOffsetCoercible = None,
                 name: str=""):

        self._j = j
        self._exchange_matrix = np.diag([0.0, j, j], dtype=float)

        super().__init__(site_1=site_1,
                         site_2=site_2,
                         cell_offset=cell_offset,
                         exchange_matrix=self._exchange_matrix,
                         name=name)

    def _exchange_serialise(self, context: SPWSerialisationContext) -> dict:
        json = self._base_serialisation(context)
        json["j"] = self._j
        return json

    @staticmethod
    def _exchange_deserialise(data: dict, context: SPWSerialisationContext):
        base = Exchange._base_deserialisation(data, context)

        return XYExchange(
            site_1=base.site_1,
            site_2=base.site_2,
            cell_offset=base.cell_offset,
            name=base.name,
            j = data["j"])

    def is_symmetric(self):
        """ Is this a symmetric exchange """
        return True

    def updated(self,
                site_1: LatticeSite | None = None,
                site_2: LatticeSite | None = None,
                cell_offset: CellOffset | None = None,
                name: str | None = None,
                j: float | None = None):
        """ Get version of this exchange with specified parameters updated"""
        return XYExchange(
                site_1=self.site_1 if site_1 is None else site_1,
                site_2=self.site_2 if site_2 is None else site_2,
                cell_offset=self.cell_offset if cell_offset is None else cell_offset,
                name=self.name if name is None else name,
                j=self.j if j is None else j)

class XXZExchange(Exchange):
    """ "XXZ" exchange, which takes the form

    H_ij = J_xy (S^x_i S^x_j + S^y_i S^y_j) + J_z (S^z_i S^z_j)

    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param j_xy: The exchange coefficient for the x and y components.
    :param j_z: The exchange coefficient for the z component

    """

    exchange_type = "XXZ"
    parameters = ["j_xy", "j_z"]
    parameter_defaults = [1.0, 1.0]
    short_string = "J"

    def __init__(self, site_1: LatticeSite, site_2: LatticeSite,
                 j_xy: float, j_z: float,
                 cell_offset: CellOffsetCoercible = None,
                 name: str=""):

        self._j_xy = j_xy
        self._j_z = j_z

        self._exchange_matrix = np.diag([j_xy, j_xy, j_z])

        super().__init__(site_1=site_1,
                         site_2=site_2,
                         cell_offset=cell_offset,
                         exchange_matrix=self._exchange_matrix,
                         name=name)



    @property
    def j_xy(self):
        """ Exchange constant for x and y """
        return self._j_xy


    @property
    def j_z(self):
        """ Exchange constant for z """
        return self._j_z

    def _exchange_serialise(self, context: SPWSerialisationContext) -> dict:
        json = self._base_serialisation(context)
        json["j_xy"] = self._j_xy
        json["j_z"] = self._j_z
        return json

    @staticmethod
    def _exchange_deserialise(data: dict, context: SPWSerialisationContext):
        base = Exchange._base_deserialisation(data, context)

        return XXZExchange(
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
        """ Get version of this exchange with specified parameters updated"""
        return XXZExchange(
                site_1=self.site_1 if site_1 is None else site_1,
                site_2=self.site_2 if site_2 is None else site_2,
                cell_offset=self.cell_offset if cell_offset is None else cell_offset,
                name=self.name if name is None else name,
                j_xy=self.j_xy if j_xy is None else j_xy,
                j_z=self.j_z if j_z is None else j_z)
    def is_symmetric(self):
        """ Is this a symmetric exchange """
        return True

class IsingExchange(Exchange):
    """Ising exchange (z component only), which takes the form

    H_ij = J_ij S^z_i S^z_j

    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param j: Scalar. The exchange coefficient J_ij.

    """

    exchange_type = "Ising"
    parameters = ["j_z"]
    parameter_defaults = [1.0]
    short_string = "J"

    def __init__(self, site_1: LatticeSite, site_2: LatticeSite,
                 j_z: float,
                 cell_offset: CellOffsetCoercible = None,
                 name: str=""):

        self._j_z = j_z

        self._exchange_matrix = np.diag([j_z])

        super().__init__(site_1=site_1,
                         site_2=site_2,
                         cell_offset=cell_offset,
                         exchange_matrix=self._exchange_matrix,
                         name=name)

    @property
    def j_z(self):
        """ Exchange constant for z """
        return self._j_z

    def _exchange_serialise(self, context: SPWSerialisationContext) -> dict:
        json = self._base_serialisation(context)
        json["j_z"] = self._j_z
        return json

    @staticmethod
    def _exchange_deserialise(data: dict, context: SPWSerialisationContext):
        base = Exchange._base_deserialisation(data, context)

        return IsingExchange(
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
        """ Get version of this exchange with specified parameters updated"""
        return IsingExchange(
                site_1=self.site_1 if site_1 is None else site_1,
                site_2=self.site_2 if site_2 is None else site_2,
                cell_offset=self.cell_offset if cell_offset is None else cell_offset,
                name=self.name if name is None else name,
                j_z=self.j_z if j_z is None else j_z)

    def is_symmetric(self):
        """ Is this a symmetric exchange """
        return True

class DMExchange(Exchange):
    """Dzyaloshinskii–Moriya exchange, which takes the form

    H_ij = D_ij (S_i x S_j)

    :param site_1: Identifier for S_i
    :param site_2: Identifier for S_j
    :param d_x: x component of the d vector above
    :param d_y: x component of the d vector above
    :param d_z: x component of the d vector above

    """

    exchange_type = "Dzyaloshinskii-Moriya"
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

        self._exchange_matrix = triple_product_matrix(np.array([d_x, d_y, d_z], dtype=float))

        super().__init__(site_1=site_1,
                         site_2=site_2,
                         cell_offset=cell_offset,
                         exchange_matrix=self._exchange_matrix,
                         name=name)


    @property
    def d_x(self):
        """ DM exchange constant for x """
        return self._d_x

    @property
    def d_y(self):
        """ DM exchange constant for y """
        return self._d_y

    @property
    def d_z(self):
        """ DM exchange constant for z """
        return self._d_z

    def _exchange_serialise(self, context: SPWSerialisationContext) -> dict:
        json = self._base_serialisation(context)
        json["d_x"] = self._d_x
        json["d_y"] = self._d_y
        json["d_z"] = self._d_z
        return json

    @staticmethod
    def _exchange_deserialise(data: dict, context: SPWSerialisationContext):
        base = Exchange._base_deserialisation(data, context)

        return DMExchange(
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
        """ Get version of this exchange with specified parameters updated"""
        return DMExchange(
                site_1=self.site_1 if site_1 is None else site_1,
                site_2=self.site_2 if site_2 is None else site_2,
                cell_offset=self.cell_offset if cell_offset is None else cell_offset,
                name=self.name if name is not None else name,
                d_x=self.d_x if d_x is None else d_x,
                d_y=self.d_y if d_y is None else d_y,
                d_z=self.d_z if d_z is None else d_z)


    def is_symmetric(self):
        """ Is this a symmetric exchange """
        return self.d_x == 0 and self.d_y == 0 and self.d_z == 0


all_exchanges = [HeisenbergExchange, DiagonalExchange, XYExchange, IsingExchange, DMExchange]
exchanges_lookup = {exchange.exchange_type: exchange for exchange in all_exchanges}
lowercase_exchange_lookup = {exchange.exchange_type.lower(): exchange for exchange in all_exchanges}
