"""Exchange terms between lattice sites."""
from dataclasses import dataclass

import numpy as np

from pyspinw.cell_offsets import CellOffsetCoercible, CellOffset
from pyspinw.checks import check_sizes
from pyspinw.exchangemetadata import ExchangeMetadata
from pyspinw.serialisation import SPWSerialisationContext, SPWSerialisable, numpy_serialise, \
    expects_keys, numpy_deserialise, SPWDeserialisationContext
from pyspinw.site import LatticeSite
from pyspinw.symmetry.group import SpaceGroup
from pyspinw.symmetry.operations import SpaceOperation
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.tolerances import tolerances
from pyspinw.util import triple_product_matrix


_exchange_id_counter = -1
def _generate_unique_exchange_id():
    """Generate a unique ID for each exchange currently loaded."""
    global _exchange_id_counter # noqa: PLW0603
    _exchange_id_counter += 1
    return _exchange_id_counter

@dataclass
class ExchangeBaseDeserialisation:
    """Hold basic properties of an exchange term."""

    name: str
    site_1: LatticeSite
    site_2: LatticeSite
    cell_offset: CellOffset
    metadata: ExchangeMetadata


class Exchange(SPWSerialisable):
    """Represent an exchange term between two lattice sites.

    Parameters
    ----------
    site_1 : LatticeSite
        First lattice site in the exchange term.
    site_2 : LatticeSite
        Second lattice site in the exchange term.
    cell_offset : CellOffsetCoercible
        Offset between the unit cells containing the two sites.
    exchange_matrix : numpy.ndarray
        3x3 matrix defining the exchange interaction.
    name : str, optional
        Name for the exchange term. Default is ``""``.
    color : tuple of float, optional
        RGB color used when displaying the exchange term.
    metadata : ExchangeMetadata, optional
        Metadata attached to the exchange term.
    """

    serialisation_name = "exchange"

    @check_sizes(exchange_matrix=(3,3))
    def __init__(self,
                 site_1: LatticeSite,
                 site_2: LatticeSite,
                 cell_offset: CellOffsetCoercible,
                 exchange_matrix: np.ndarray,
                 name: str = "",
                 color: tuple[float, float, float] | None = None,
                 metadata: ExchangeMetadata | None = None):

        self._name = name
        self._site_1 = site_1
        self._site_2 = site_2
        self._cell_offset = CellOffset.coerce(cell_offset)
        self._exchange_matrix = exchange_matrix

        self._unique_id = _generate_unique_exchange_id()
        self.metadata = ExchangeMetadata() if metadata is None else metadata

        if color is not None:
            self.metadata.color = color


    exchange_type = "General"
    parameters: list[str] = []
    parameter_defaults: list[float] = []
    short_string = "M"


    @property
    def name(self):
        """Name of this exchange."""
        return self._name

    @property
    def unique_id(self) -> int:
        """Unique identifier for this exchange."""
        return self._unique_id

    @property
    def site_1(self):
        """First site in the exchange term."""
        return self._site_1

    @property
    def site_2(self):
        """Second site in the exchange term."""
        return self._site_2

    @property
    def cell_offset(self):
        """Offset between unit cells containing sites 1 and 2."""
        return self._cell_offset

    @property
    def exchange_matrix(self) -> np.ndarray:
        r"""The exchange matrix for this exchange.

        If :math:`H` is the energy contribution for this exchange,
        :math:`\mathbf{S}_i` and :math:`\mathbf{S}_j` are the spin states,
        and :math:`M` is the exchange matrix, then

        .. math::

            H = \mathbf{S}_i^T M \mathbf{S}_j.
        """
        return self._exchange_matrix

    @property
    def parameter_string(self) -> str:
        """String representation of the exchange parameters."""
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
        """Vector from site 1 to site 2 in lattice coordinates."""
        return self.cell_offset.vector + self._site_2.ijk - self._site_1.ijk

    def vector(self, unit_cell: UnitCell):
        """Return the vector from site 1 to site 2 in Cartesian coordinates.

        Parameters
        ----------
        unit_cell
            Unit cell used to convert lattice coordinates to Cartesian coordinates.

        Returns
        -------
        numpy.ndarray
            Vector from site 1 to site 2 in Cartesian coordinates.
        """
        return unit_cell.lattice_units_to_cartesian(self.lattice_vector)

    def distance(self, unit_cell: UnitCell):
        """Return the distance between site 1 and site 2 in Cartesian coordinates.

        Parameters
        ----------
        unit_cell
            Unit cell used to convert lattice coordinates to Cartesian coordinates.

        Returns
        -------
        float
            Distance between site 1 and site 2 in Cartesian coordinates.
        """
        return np.sqrt(np.sum(self.vector(unit_cell)**2))

    def _base_serialisation(self, context: SPWSerialisationContext):
        return {
            "name": self.name,
            "site_1": self._site_1._serialise(context),
            "site_2": self._site_2._serialise(context),
            "cell_offset": self._cell_offset._serialise(context),
            "metadata": self.metadata._serialise(context)
        }

    @check_sizes(exchange_matrix=(3,3), force_numpy=True, allow_nones=True)
    def updated(self,
                site_1 :LatticeSite | None = None,
                site_2: LatticeSite | None = None,
                cell_offset: CellOffset | None = None,
                name: str | None = None,
                exchange_matrix: np.ndarray | None = None,
                metadata: ExchangeMetadata | None = None):
        """Return a copy of this exchange with variables replaced.

        Parameters
        ----------
        site_1 : LatticeSite, optional
            Replacement first lattice site. If omitted, the current site is reused.
        site_2 : LatticeSite, optional
            Replacement second lattice site. If omitted, the current site is reused.
        cell_offset : CellOffset, optional
            Replacement unit-cell offset. If omitted, the current offset is reused.
        name : str, optional
            Replacement exchange name. If omitted, the current name is reused.
        exchange_matrix : numpy.ndarray, optional
            Replacement exchange matrix. If omitted, the current matrix is reused.
        metadata : ExchangeMetadata, optional
            Replacement metadata. If omitted, the current metadata is copied.
        """
        return Exchange(
            site_1=self.site_1 if site_1 is None else site_1,
            site_2=self.site_2 if site_2 is None else site_2,
            cell_offset=self.cell_offset if cell_offset is None else cell_offset,
            name=self.name if name is None else name,
            exchange_matrix=self.exchange_matrix if exchange_matrix is None else exchange_matrix,
            metadata=self.metadata.copy() if metadata is None else metadata.copy())

    @staticmethod
    @expects_keys("name, site_1, site_2, cell_offset")
    def _base_deserialisation(json, context: SPWSerialisationContext):
        return ExchangeBaseDeserialisation(
            name=json["name"],
            site_1=LatticeSite._deserialise(json["site_1"], context),
            site_2=LatticeSite._deserialise(json["site_2"], context),
            cell_offset=CellOffset._deserialise(json["cell_offset"], context),
            metadata=ExchangeMetadata._deserialise(json["metadata"], context)
        )

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
        """Return whether this is a symmetric exchange."""
        return np.all(np.abs(self.exchange_matrix - self.exchange_matrix.T) < tolerances.IS_ZERO_TOL)

    def _obeys_symmetry(self,
                        unit_cell: UnitCell,
                        identity_operations: set[SpaceOperation],
                        inversion_operations: set[SpaceOperation]) -> bool:

        """ Main logic for symmetry checking """

        # We want the exchange matrix in lattice units TODO: Verify the details of this transform
        exchange_matrix = unit_cell._xyz_spins @ self._exchange_matrix @ unit_cell._xyz_spins.T

        for operation in identity_operations:
            if not np.allclose(exchange_matrix,
                           operation.point_operation_matrix @ exchange_matrix @ operation.point_operation_matrix.T):
                return False

        exchange_matrix_T = self._exchange_matrix.T
        for operation in inversion_operations:
            if not np.allclose(exchange_matrix,
                               operation.point_operation_matrix @ exchange_matrix_T @
                               operation.point_operation_matrix.T):
                return False

        return True

    def obeys_symmetry(self, structure: "Structure") -> bool:
        """ Check that this exchange is consistent with the symmetry group """
        spacegroup = structure.spacegroup
        unit_cell = structure.unit_cell

        # Checking is easier than finding the list of symmetry groups
        identity_operations, inversion_operations = spacegroup.operations_on_single_site_pairs(self.site_1, self.site_2)
        return self._obeys_symmetry(unit_cell, identity_operations, inversion_operations)

    def symmetry_copy(self,
                      structure: "Structure",
                      site_1: LatticeSite,
                      site_2: LatticeSite,
                      cell_offset: CellOffsetCoercible = (0,0,0)):
        """ Copy this exchange using symmetry operations """
        # We want to copy the exchange under symmetry operations
        # There might be more than one symmetry operation that maps the pair of sites
        #  however, the effect on the exchange should be the same for all these operations,
        #  this means we can just pick an arbitrary one.
        # If this turns out not to be the case, then exchange itself does not need to obey the
        # symmetry constraints

        spacegroup = structure.spacegroup
        unit_cell = structure.unit_cell

        if self.obeys_symmetry(spacegroup):
            # find the operations that map the pairs

            pair_operations = spacegroup.operations_between_pairs(
                (self.site_1, self.site_2),
                (site_1, site_2))

            if len(pair_operations) == 0:
                raise ValueError("New points are not related to the original by symmetry")

            # Pick one element for the transformation
            op = next(iter(pair_operations))
            transform = op.point_operation_matrix

            # Apply after transforming to xyz space TODO: Verify, could require inverses/transforms
            exchange_matrix_ijk = unit_cell._xyz_spins @ self.exchange_matrix @ unit_cell._xyz_spins.T
            new_exchange_matrix_ijk = transform @ exchange_matrix_ijk @ transform.T
            new_exchange_matrix = unit_cell._xyz_spins_inv @ new_exchange_matrix_ijk @ unit_cell._xyz_spins_inv.T

            return Exchange(site_1, site_2,
                            exchange_matrix=new_exchange_matrix,
                            name=f"{self.name} [{op.text_form}]",
                            cell_offset=CellOffset.coerce(cell_offset))

        else:
            raise ValueError("Exchange does not obey symmetry constraints, cannot use symmetry to copy")

    def symmetry_fill(self, structure: "Structure"):
        """ Make multiple copies of this exchange so that symmetry is satisfied """

class HeisenbergExchange(Exchange):
    r"""Represent a Heisenberg exchange term.

    The exchange takes the form

    .. math::

        H_{ij} = J_{ij} \, (\mathbf{S}_i \cdot \mathbf{S}_j).

    Parameters
    ----------
    site_1 : LatticeSite
        Lattice site associated with :math:`\mathbf{S}_i`.
    site_2 : LatticeSite
        Lattice site associated with :math:`\mathbf{S}_j`.
    j : float
        Exchange coefficient, :math:`J_{ij}`.
    cell_offset : CellOffsetCoercible, optional
        Offset between the unit cells containing the two sites.
    name : str, optional
        Name for the exchange term. Default is ``""``.
    color : tuple of float, optional
        RGB color used when displaying the exchange term.
    metadata : ExchangeMetadata, optional
        Metadata attached to the exchange term.
    """

    exchange_type = "Heisenberg"
    parameters = ["j"]
    parameter_defaults = [1.0]
    short_string = "J"


    def __init__(self, site_1: LatticeSite, site_2: LatticeSite,
                 j: float,
                 cell_offset: CellOffsetCoercible = None,
                 name: str="",
                 color: tuple[float, float, float] | None = None,
                 metadata: ExchangeMetadata | None = None):

        self._j = j
        self._exchange_matrix = j * np.eye(3)

        super().__init__(site_1=site_1,
                         site_2=site_2,
                         cell_offset=cell_offset,
                         exchange_matrix=self._exchange_matrix,
                         name=name,
                         color=color,
                         metadata=metadata)

    @property
    def j(self):
        """Exchange constant."""
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
        """Return whether this is a symmetric exchange.

        A Heisenberg exchange is always symmetric, so this always returns ``True``.
        """
        return True

    def updated(self,
                site_1 :LatticeSite | None = None,
                site_2: LatticeSite | None = None,
                cell_offset: CellOffset | None = None,
                name: str | None = None,
                j: float | None = None,
                metadata: ExchangeMetadata | None = None
                ):
        """Return a copy of this exchange with variables replaced.

        Parameters
        ----------
        site_1 : LatticeSite, optional
            Replacement first lattice site. If omitted, the current site is reused.
        site_2 : LatticeSite, optional
            Replacement second lattice site. If omitted, the current site is reused.
        cell_offset : CellOffset, optional
            Replacement unit-cell offset. If omitted, the current offset is reused.
        name : str, optional
            Replacement exchange name. If omitted, the current name is reused.
        j : float, optional
            Replacement exchange coefficient. If omitted, the current coefficient is reused.
        metadata : ExchangeMetadata, optional
            Replacement metadata. If omitted, the current metadata is copied.
        """
        return HeisenbergExchange(
                site_1=self.site_1 if site_1 is None else site_1,
                site_2 = self.site_2 if site_2 is None else site_2,
                cell_offset = self.cell_offset if cell_offset is None else cell_offset,
                name = self.name if name is None else name,
                j = self.j if j is None else j,
                metadata=self.metadata.copy() if metadata is None else metadata.copy()
                )




class DiagonalExchange(Exchange):
    r"""Represent a diagonal exchange term.

    The exchange takes the form

    .. math::

        H_{ij} = J^{xx}_{ij} \mathbf{S}^x_i \mathbf{S}^x_j + J^{yy}_{ij} \mathbf{S}^y_i \mathbf{S}^y_j
                 + J^{zz}_{ij} \mathbf{S}^z_i \mathbf{S}^z_j.

    Parameters
    ----------
    site_1 : LatticeSite
        Lattice site associated with :math:`\mathbf{S}_i`.
    site_2 : LatticeSite
        Lattice site associated with :math:`\mathbf{S}_j`.
    j_x : float
        Scalar exchange coefficient for the x component, :math:`J^{xx}_{ij}`.
    j_y : float
        Scalar exchange coefficient for the y component, :math:`J^{yy}_{ij}`.
    j_z : float
        Scalar exchange coefficient for the z component, :math:`J^{zz}_{ij}`.
    cell_offset : CellOffsetCoercible, optional
        Offset between the unit cells containing the two sites.
    name : str, optional
        Name for the exchange term. Default is ``""``.
    color : tuple of float, optional
        RGB color used when displaying the exchange term.
    metadata : ExchangeMetadata, optional
        Metadata attached to the exchange term.
    """

    parameters = ["j_x", "j_y", "j_z"]
    parameter_defaults = [1.0, 1.0, 1.0]
    short_string = "J"

    def __init__(self, site_1: LatticeSite, site_2: LatticeSite,
                 j_x: float, j_y: float, j_z: float,
                 cell_offset: CellOffsetCoercible = None,
                 name: str="",
                 color: tuple[float, float, float] | None = None,
                 metadata: ExchangeMetadata | None = None):

        self._j_x = j_x
        self._j_y = j_y
        self._j_z = j_z

        self._exchange_matrix = np.diag([j_x, j_y, j_z])

        super().__init__(site_1=site_1,
                         site_2=site_2,
                         cell_offset=cell_offset,
                         exchange_matrix=self._exchange_matrix,
                         name=name,
                         color=color,
                         metadata=metadata)

    @property
    def j_x(self):
        """Exchange constant for x."""
        return self._j_x


    @property
    def j_y(self):
        """Exchange constant for y."""
        return self._j_y


    @property
    def j_z(self):
        """Exchange constant for z."""
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
        """Return whether this is a symmetric exchange.

        A diagonal exchange is always symmetric, so this always returns ``True``.
        """
        return True

    def updated(self,
                site_1 :LatticeSite | None = None,
                site_2: LatticeSite | None = None,
                cell_offset: CellOffset | None = None,
                name: str | None = None,
                j_x: float | None = None,
                j_y: float | None = None,
                j_z: float | None = None,
                metadata: ExchangeMetadata | None = None,
                ):
        """Return a copy of this exchange with variables replaced.

        Parameters
        ----------
        site_1 : LatticeSite, optional
            Replacement first lattice site. If omitted, the current site is reused.
        site_2 : LatticeSite, optional
            Replacement second lattice site. If omitted, the current site is reused.
        cell_offset : CellOffset, optional
            Replacement unit-cell offset. If omitted, the current offset is reused.
        name : str, optional
            Replacement exchange name. If omitted, the current name is reused.
        j_x : float, optional
            Replacement x exchange coefficient. If omitted, the current coefficient is reused.
        j_y : float, optional
            Replacement y exchange coefficient. If omitted, the current coefficient is reused.
        j_z : float, optional
            Replacement z exchange coefficient. If omitted, the current coefficient is reused.
        metadata : ExchangeMetadata, optional
            Replacement metadata. If omitted, the current metadata is copied.
        """
        return DiagonalExchange(
                site_1=self.site_1 if site_1 is None else site_1,
                site_2 = self.site_2 if site_2 is None else site_2,
                cell_offset = self.cell_offset if cell_offset is None else cell_offset,
                name = self.name if name is None else name,
                j_x = self.j_x if j_x is None else j_x,
                j_y = self.j_y if j_y is None else j_y,
                j_z = self.j_z if j_z is None else j_z,
                metadata=self.metadata.copy() if metadata is None else metadata.copy())

class XYExchange(Exchange):
    r"""Represent an XY exchange term.

    The exchange takes the form

    .. math::

        H_{ij} = J_{ij} \, (\mathbf{S}^x_i \mathbf{S}^x_j + \mathbf{S}^y_i \mathbf{S}^y_j).

    Parameters
    ----------
    site_1 : LatticeSite
        Lattice site associated with :math:`\mathbf{S}_i`.
    site_2 : LatticeSite
        Lattice site associated with :math:`\mathbf{S}_j`.
    j : float
        Exchange coefficient for the x and y components, :math:`J_{ij}`.
    cell_offset : CellOffsetCoercible, optional
        Offset between the unit cells containing the two sites.
    name : str, optional
        Name for the exchange term. Default is ``""``.
    color : tuple of float, optional
        RGB color used when displaying the exchange term.
    metadata : ExchangeMetadata, optional
        Metadata attached to the exchange term.
    """

    exchange_type = "XY"
    parameters = ["j"]
    parameter_defaults = [1.0]
    short_string = "J"


    @property
    def j(self):
        """Exchange constant for x and y."""
        return self._j

    def __init__(self, site_1: LatticeSite, site_2: LatticeSite,
                 j: float,
                 cell_offset: CellOffsetCoercible = None,
                 name: str="",
                 color: tuple[float, float, float] | None = None,
                 metadata: ExchangeMetadata | None = None):

        self._j = j
        self._exchange_matrix = np.diag([j, j, 0.0], dtype=float)

        super().__init__(site_1=site_1,
                         site_2=site_2,
                         cell_offset=cell_offset,
                         exchange_matrix=self._exchange_matrix,
                         name=name,
                         color=color,
                         metadata=metadata)

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
        """Return whether this is a symmetric exchange.

        An XY exchange is always symmetric, so this always returns ``True``.
        """
        return True

    def updated(self,
                site_1: LatticeSite | None = None,
                site_2: LatticeSite | None = None,
                cell_offset: CellOffset | None = None,
                name: str | None = None,
                j: float | None = None,
                metadata: ExchangeMetadata | None = None):
        """Return a copy of this exchange with variables replaced.

        Parameters
        ----------
        site_1 : LatticeSite, optional
            Replacement first lattice site. If omitted, the current site is reused.
        site_2 : LatticeSite, optional
            Replacement second lattice site. If omitted, the current site is reused.
        cell_offset : CellOffset, optional
            Replacement unit-cell offset. If omitted, the current offset is reused.
        name : str, optional
            Replacement exchange name. If omitted, the current name is reused.
        j : float, optional
            Replacement exchange coefficient. If omitted, the current coefficient is reused.
        metadata : ExchangeMetadata, optional
            Replacement metadata. If omitted, the current metadata is copied.
        """
        return XYExchange(
                site_1=self.site_1 if site_1 is None else site_1,
                site_2=self.site_2 if site_2 is None else site_2,
                cell_offset=self.cell_offset if cell_offset is None else cell_offset,
                name=self.name if name is None else name,
                j=self.j if j is None else j,
                metadata=self.metadata.copy() if metadata is None else metadata.copy())

class XXZExchange(Exchange):
    r"""Represent an XXZ exchange term.

    The exchange takes the form

    .. math::

        H_{ij} = J_{xy} \, (\mathbf{S}^x_i \mathbf{S}^x_j
        + \mathbf{S}^y_i \mathbf{S}^y_j) + J_z \, (\mathbf{S}^z_i \mathbf{S}^z_j).

    Parameters
    ----------
    site_1 : LatticeSite
        Lattice site associated with :math:`\mathbf{S}_i`.
    site_2 : LatticeSite
        Lattice site associated with :math:`\mathbf{S}_j`.
    j_xy : float
        Exchange coefficient for the x and y components, :math:`J_{xy}`.
    j_z : float
        Exchange coefficient for the z component, :math:`J_z`.
    cell_offset : CellOffsetCoercible, optional
        Offset between the unit cells containing the two sites.
    name : str, optional
        Name for the exchange term. Default is ``""``.
    color : tuple of float, optional
        RGB color used when displaying the exchange term.
    metadata : ExchangeMetadata, optional
        Metadata attached to the exchange term.
    """

    exchange_type = "XXZ"
    parameters = ["j_xy", "j_z"]
    parameter_defaults = [1.0, 1.0]
    short_string = "J"

    def __init__(self, site_1: LatticeSite, site_2: LatticeSite,
                 j_xy: float, j_z: float,
                 cell_offset: CellOffsetCoercible = None,
                 name: str="",
                 color: tuple[float, float, float] | None = None,
                 metadata: ExchangeMetadata | None = None):

        self._j_xy = j_xy
        self._j_z = j_z

        self._exchange_matrix = np.diag([j_xy, j_xy, j_z])

        super().__init__(site_1=site_1,
                         site_2=site_2,
                         cell_offset=cell_offset,
                         exchange_matrix=self._exchange_matrix,
                         name=name,
                         color=color,
                         metadata=metadata)



    @property
    def j_xy(self):
        """Exchange constant for x and y."""
        return self._j_xy


    @property
    def j_z(self):
        """Exchange constant for z."""
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
                metadata: ExchangeMetadata | None = None
                ):
        """Return a copy of this exchange with variables replaced.

        Parameters
        ----------
        site_1 : LatticeSite, optional
            Replacement first lattice site. If omitted, the current site is reused.
        site_2 : LatticeSite, optional
            Replacement second lattice site. If omitted, the current site is reused.
        cell_offset : CellOffset, optional
            Replacement unit-cell offset. If omitted, the current offset is reused.
        name : str, optional
            Replacement exchange name. If omitted, the current name is reused.
        j_xy : float, optional
            Replacement x-y exchange coefficient. If omitted, the current coefficient is reused.
        j_z : float, optional
            Replacement z exchange coefficient. If omitted, the current coefficient is reused.
        metadata : ExchangeMetadata, optional
            Replacement metadata. If omitted, the current metadata is copied.
        """
        return XXZExchange(
                site_1=self.site_1 if site_1 is None else site_1,
                site_2=self.site_2 if site_2 is None else site_2,
                cell_offset=self.cell_offset if cell_offset is None else cell_offset,
                name=self.name if name is None else name,
                j_xy=self.j_xy if j_xy is None else j_xy,
                j_z=self.j_z if j_z is None else j_z,
                metadata=self.metadata.copy() if metadata is None else metadata.copy())
    def is_symmetric(self):
        """Return whether this is a symmetric exchange.

        An XXZ exchange is always symmetric, so this always returns ``True``.
        """
        return True

class IsingExchange(Exchange):
    r"""Represent an Ising exchange term for the z component.

    The exchange takes the form

    .. math::

        H_{ij} = J_{ij} \, \mathbf{S}^z_i \mathbf{S}^z_j.

    Parameters
    ----------
    site_1 : LatticeSite
        Lattice site associated with :math:`\mathbf{S}_i`.
    site_2 : LatticeSite
        Lattice site associated with :math:`\mathbf{S}_j`.
    j_z : float
        Exchange coefficient for the z component, :math:`J_{ij}`.
    cell_offset : CellOffsetCoercible, optional
        Offset between the unit cells containing the two sites.
    name : str, optional
        Name for the exchange term. Default is ``""``.
    color : tuple of float, optional
        RGB color used when displaying the exchange term.
    metadata : ExchangeMetadata, optional
        Metadata attached to the exchange term.
    """

    exchange_type = "Ising"
    parameters = ["j_z"]
    parameter_defaults = [1.0]
    short_string = "J"

    def __init__(self, site_1: LatticeSite, site_2: LatticeSite,
                 j_z: float,
                 cell_offset: CellOffsetCoercible = None,
                 name: str="",
                 color: tuple[float, float, float] | None = None,
                 metadata: ExchangeMetadata | None = None):
        self._j_z = j_z

        self._exchange_matrix = np.diag([0.0, 0.0, j_z])

        super().__init__(site_1=site_1,
                         site_2=site_2,
                         cell_offset=cell_offset,
                         exchange_matrix=self._exchange_matrix,
                         name=name,
                         color=color,
                         metadata=metadata)

    @property
    def j_z(self):
        """Exchange constant for z."""
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
                j_z: float | None = None,
                metadata: ExchangeMetadata | None = None):
        """Return a copy of this exchange with variables replaced.

        Parameters
        ----------
        site_1 : LatticeSite, optional
            Replacement first lattice site. If omitted, the current site is reused.
        site_2 : LatticeSite, optional
            Replacement second lattice site. If omitted, the current site is reused.
        cell_offset : CellOffset, optional
            Replacement unit-cell offset. If omitted, the current offset is reused.
        name : str, optional
            Replacement exchange name. If omitted, the current name is reused.
        j_z : float, optional
            Replacement z exchange coefficient. If omitted, the current coefficient is reused.
        metadata : ExchangeMetadata, optional
            Replacement metadata. If omitted, the current metadata is copied.
        """
        return IsingExchange(
                site_1=self.site_1 if site_1 is None else site_1,
                site_2=self.site_2 if site_2 is None else site_2,
                cell_offset=self.cell_offset if cell_offset is None else cell_offset,
                name=self.name if name is None else name,
                j_z=self.j_z if j_z is None else j_z,
                metadata=self.metadata.copy() if metadata is None else metadata.copy())

    def is_symmetric(self):
        """Return whether this is a symmetric exchange.

        An Ising exchange is always symmetric, so this always returns ``True``.
        """
        return True

class DMExchange(Exchange):
    r"""Represent a Dzyaloshinskii-Moriya exchange term.

    The exchange takes the form

    .. math::

        H_{ij} = \mathbf{D}_{ij} \cdot (\mathbf{S}_i \times \mathbf{S}_j).

    Parameters
    ----------
    site_1 : LatticeSite
        Lattice site associated with :math:`\mathbf{S}_i`.
    site_2 : LatticeSite
        Lattice site associated with :math:`\mathbf{S}_j`.
    d_x : float
        x component of the D vector, :math:`\mathbf{D}^x_{ij}`.
    d_y : float
        y component of the D vector, :math:`\mathbf{D}^y_{ij}`.
    d_z : float
        z component of the D vector, :math:`\mathbf{D}^z_{ij}`.
    cell_offset : CellOffsetCoercible, optional
        Offset between the unit cells containing the two sites.
    name : str, optional
        Name for the exchange term. Default is ``""``.
    color : tuple of float, optional
        RGB color used when displaying the exchange term.
    metadata : ExchangeMetadata, optional
        Metadata attached to the exchange term.
    """

    exchange_type = "Dzyaloshinskii-Moriya"
    parameters = ["d_x", "d_y", "d_z"]
    parameter_defaults = [1.0, 1.0, 1.0]
    short_string = "DM"


    def __init__(self, site_1: LatticeSite, site_2: LatticeSite,
                 d_x: float, d_y: float, d_z: float,
                 cell_offset: CellOffsetCoercible = None,
                 name: str="",
                 color: tuple[float, float, float] | None = None,
                 metadata: ExchangeMetadata | None = None):

        self._d_x = d_x
        self._d_y = d_y
        self._d_z = d_z

        self._exchange_matrix = triple_product_matrix(np.array([d_x, d_y, d_z], dtype=float))

        super().__init__(site_1=site_1,
                         site_2=site_2,
                         cell_offset=cell_offset,
                         exchange_matrix=self._exchange_matrix,
                         name=name,
                         color=color,
                         metadata=metadata)


    @property
    def d_x(self):
        """DM exchange constant for x."""
        return self._d_x

    @property
    def d_y(self):
        """DM exchange constant for y."""
        return self._d_y

    @property
    def d_z(self):
        """DM exchange constant for z."""
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
                d_z: float | None = None,
                metadata: ExchangeMetadata | None = None):
        """Return a copy of this exchange with variables replaced.

        Parameters
        ----------
        site_1 : LatticeSite, optional
            Replacement first lattice site. If omitted, the current site is reused.
        site_2 : LatticeSite, optional
            Replacement second lattice site. If omitted, the current site is reused.
        cell_offset : CellOffset, optional
            Replacement unit-cell offset. If omitted, the current offset is reused.
        name : str, optional
            Replacement exchange name. If omitted, the current name is reused.
        d_x : float, optional
            Replacement x D-vector component. If omitted, the current component is reused.
        d_y : float, optional
            Replacement y D-vector component. If omitted, the current component is reused.
        d_z : float, optional
            Replacement z D-vector component. If omitted, the current component is reused.
        metadata : ExchangeMetadata, optional
            Replacement metadata. If omitted, the current metadata is copied.
        """
        return DMExchange(
                site_1=self.site_1 if site_1 is None else site_1,
                site_2=self.site_2 if site_2 is None else site_2,
                cell_offset=self.cell_offset if cell_offset is None else cell_offset,
                name=self.name if name is not None else name,
                d_x=self.d_x if d_x is None else d_x,
                d_y=self.d_y if d_y is None else d_y,
                d_z=self.d_z if d_z is None else d_z,
                metadata=self.metadata.copy() if metadata is None else metadata.copy())


    def is_symmetric(self):
        """Return whether this is a symmetric exchange."""
        return self.d_x == 0 and self.d_y == 0 and self.d_z == 0


all_exchanges = [HeisenbergExchange, DiagonalExchange, XYExchange, IsingExchange, DMExchange]
exchanges_lookup = {exchange.exchange_type: exchange for exchange in all_exchanges}
lowercase_exchange_lookup = {exchange.exchange_type.lower(): exchange for exchange in all_exchanges}
