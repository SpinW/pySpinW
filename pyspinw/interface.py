"""Helper functions for python interface."""
import math

import numpy as np
from numpy._typing import ArrayLike

from pyspinw.anisotropy import AxisMagnitudeAnisotropy, Anisotropy
from pyspinw.batch_exchanges import default_naming_pattern
from pyspinw.checks import check_sizes
from pyspinw.exchange import Exchange, HeisenbergExchange
from pyspinw.exchangegroup import DirectionalityFilter, InPlaneFilter, InDirectionFilter, ExchangeGroup, \
    BiDirectionFilter
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.site import Site
from pyspinw.structure import Structure
from pyspinw.symmetry.group import database, NoSuchGroup, ExactMatch, PartialMatch
from pyspinw.symmetry.supercell import PropagationVector, CommensuratePropagationVector, RotationTransform, \
    TransformationSupercell, SummationSupercell, RotationSupercell
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.units import CoordsUnits
from pyspinw.tb2j import TB2J_Input

try:
    from pyspinw.gui.viewer import show_object
except ModuleNotFoundError:
    def show_hamiltonian(*args, **kwargs):
        """Show a Hamiltonian in the viewer."""
        raise RuntimeError('PySide or OpenGL not installed')

def _check_positions_spins_shape(positions: ArrayLike,
                                   spins: ArrayLike | None = None) -> tuple[ArrayLike, ArrayLike]:
    """Check position and spin arrays and put them into canonical shape.

    Parameters
    ----------
    positions : array_like
        Atomic site positions as a length-3 vector or ``(n_sites, 3)`` array.
    spins : array_like, optional
        Site spin vectors as a length-3 vector, ``(n_sites, 3)`` array, or
        ``(n_sites, m, 3)`` array. If ``None``, all spins are initialised to zero.
        Default is ``None``.

    Returns
    -------
    positions : (n_sites, 3) numpy.ndarray
        Canonical 2D array of site positions.
    spins : numpy.ndarray
        Canonical 3D array of site spin components.

    Raises
    ------
    ValueError
        If ``positions`` or ``spins`` do not have valid dimensions or matching
        site counts.
    """
    # check the positions are n-by-3, or 1d
    positions = np.array(positions, dtype=float)
    match len(positions.shape):
        case 1:
            if positions.shape[0] == 3:
                positions = positions.reshape(1, 3)
                n_sites = 1
            else:
                raise ValueError("Expected positions to be n-by-3, or a single vector of length 3")
        case 2:
            if positions.shape[1] == 3:
                n_sites = positions.shape[0]
            else:
                raise ValueError("Expected positions to be n-by-3, or a single vector of length 3")
        case _:
            raise ValueError("Expected positions to be n-by-3, or a single vector of length 3")

    # check spin definitions
    if spins is None:
        spins = np.zeros((n_sites, 3, 1))
    else:
        spins = np.array(spins, dtype=float)
    match len(spins.shape):
        case 1:
            if n_sites != 1:
                raise ValueError("Expected same number of spins as positions")

            if spins.shape[0] != 3:
                raise ValueError("Expected spin to have three components")

            spins = spins.reshape(1,3,1)

        case 2:
            if n_sites != spins.shape[0]:
                raise ValueError("Expected same number of spins as positions")

            if spins.shape[1] != 3:
                raise ValueError("Expected spins to have three components")

            spins = spins.reshape(-1, 1, 3)

        case 3:
            if n_sites != spins.shape[0]:
                raise ValueError("Expected same number of spins as postions")

            if spins.shape[2] != 3:
                raise ValueError("Expected spins to have three components")

        case _:
            raise ValueError("Expected spins to be n-by-3-m, n-by-3, or vector of length 3")

    return positions, spins

def generate_sites(positions: ArrayLike,
                   spins: ArrayLike | None = None,
                   names: list[str] | None=None) -> list[Site]:
    """Create a list of lattice sites.

    Parameters
    ----------
    positions : array_like
        Positions of the sites in fractional lattice coordinates.
    spins : array_like, optional
        Spins of the sites; if not specified, they will be set to zero.
        Default is ``None``.
    names : list of str, optional
        A list of names for sites. Default is ``None``.

    Returns
    -------
    list of Site
        Created ``Site`` objects.
    """
    positions, spins = _check_positions_spins_shape(positions, spins)

    n_sites = positions.shape[0]
    if names is None:
        names = [f'site{i:02d}' for i in range(n_sites)]

    out = []
    for i in range(n_sites):
        position = positions[i, :]

        out.append(Site(
            i = position[0],
            j = position[1],
            k = position[2],
            supercell_spins=np.squeeze(spins[i, :, :]),
            name = names[i]
        ))

    return out

def _transform_site(unit_cell: UnitCell,
                    positions: ArrayLike,
                    spins: ArrayLike,
                    names: list[str] | None,
                    magnitudes: ArrayLike | None,
                    positions_unit: CoordsUnits | str,
                    spins_unit: CoordsUnits | str) -> tuple:
    """Convert site data into default units.

    Parameters
    ----------
    unit_cell : UnitCell
        Unit cell used for coordinate conversions.
    positions : array_like
        Positions of the sites.
    spins : array_like
        Spin vectors of the sites.
    names : list of str or None
        Names for the sites.
    magnitudes : array_like or None
        Target spin magnitudes for each site.
    positions_unit : CoordsUnits or str
        Coordinate unit of ``positions`` (``'lu'`` or ``'xyz'``).
    spins_unit : CoordsUnits or str
        Coordinate unit of ``spins`` (``'lu'`` or ``'xyz'``).

    Returns
    -------
    positions : numpy.ndarray
        Site positions in fractional lattice units.
    spins : numpy.ndarray
        Site spin vectors in Cartesian coordinates scaled to ``magnitudes``.
    names : list of str or None
        Names for the sites.
    """
    positions, spins = _check_positions_spins_shape(positions, spins)
    # check units
    positions_unit = CoordsUnits(positions_unit)
    spins_unit = CoordsUnits(spins_unit)
    # Convert cell positions
    if positions_unit == CoordsUnits.XYZ:
        positions = unit_cell.cartesian_to_lattice_units(np.array(positions))
    # convert spins
    if spins_unit == CoordsUnits.LU:
        spins = np.array(spins, dtype=float)
        for i in range(spins.shape[0]):
            magnitude = np.linalg.norm(spins[i,:,:]) if magnitudes is None else magnitudes[i]
            spin = unit_cell.lattice_units_to_cartesian(spins[i, :, :])
            spins[i,:,:] = spin * (magnitude / np.linalg.norm(spin))
    elif magnitudes is not None:
        for i in range(spins.shape[0]):
            spins[i,:,:] *= (magnitudes[i] / np.linalg.norm(spins[i,:,:]))
    return positions, spins, names

def generate_structure(unit_cell: UnitCell,
                       positions: ArrayLike,
                       spins: ArrayLike,
                       propagation_vectors: ArrayLike | None = None,
                       names: list[str] | None = None,
                       magnitudes: ArrayLike | None = None,
                       positions_unit: CoordsUnits | str = 'lu',
                       spins_unit: CoordsUnits | str = 'xyz') -> Structure:
    """Create a magnetic structure from a list of positions and spins.

    Parameters
    ----------
    unit_cell : UnitCell
        A ``UnitCell`` object.
    positions : array_like
        Positions of the sites.
    spins : array_like
        Spins of the sites.
    propagation_vectors : array_like, optional
        The propagation vector(s). Default is ``None``.
    names : list of str, optional
        A list of names for sites. Default is ``None``.
    magnitudes : array_like, optional
        A list of spin magnitudes (spin length ``S``). If not specified, the norm
        of the spin vectors will be used as the spin length. Default is ``None``.
    positions_unit : CoordsUnits or str, optional
        The units for atomic positions; either ``'lu'`` or ``'xyz'``. Default is
        ``'lu'`` (lattice units).
    spins_unit : CoordsUnits or str, optional
        The units for the magnetic spin directions; either ``'lu'`` or ``'xyz'``.
        Default is ``'xyz'`` (Cartesian axis with ``x || a``).

    Returns
    -------
    Structure
        Constructed magnetic ``Structure``.
    """
    s = generate_sites(*_transform_site(unit_cell, positions, spins, names, magnitudes, positions_unit, spins_unit))
    if propagation_vectors is not None:
        pvs = np.array(propagation_vectors, ndmin=2)
        if pvs.shape[0] != 3 and pvs.shape[1] == 3:
            pvs = pvs.T
        pvs = [CommensuratePropagationVector(*pvs[:,ip]) for ip in range(pvs.shape[1])]
        return Structure(s, unit_cell, supercell=SummationSupercell(propagation_vectors=pvs))
    return Structure(s, unit_cell)

def generate_helical_structure(unit_cell: UnitCell,
                               positions: ArrayLike,
                               spins: ArrayLike,
                               perpendicular: ArrayLike,
                               propagation_vector: ArrayLike,
                               names: list[str] | None=None,
                               magnitudes: ArrayLike | None = None,
                               positions_unit: CoordsUnits | str = 'lu',
                               spins_unit: CoordsUnits | str = 'xyz') -> Structure:
    """Create a helical structure with a propagation vector and plane normal.

    Parameters
    ----------
    unit_cell : UnitCell
        A ``UnitCell`` object.
    positions : array_like
        Positions of the sites.
    spins : array_like
        Spins of the sites.
    perpendicular : array_like
        Vector perpendicular to the helix plane of rotation.
    propagation_vector : array_like
        The propagation vector.
    names : list of str, optional
        A list of names for sites. Default is ``None``.
    magnitudes : array_like, optional
        A list of spin magnitudes (spin length ``S``). If not specified, the norm
        of the spin vectors will be used as the spin length. Default is ``None``.
    positions_unit : CoordsUnits or str, optional
        The units for atomic positions; either ``'lu'`` or ``'xyz'``. Default is
        ``'lu'`` (lattice units).
    spins_unit : CoordsUnits or str, optional
        The units for the magnetic spin directions; either ``'lu'`` or ``'xyz'``.
        Default is ``'xyz'`` (Cartesian axis with ``x || a``).

    Returns
    -------
    Structure
        Constructed helical ``Structure``.
    """
    s = generate_sites(*_transform_site(unit_cell, positions, spins, names, magnitudes, positions_unit, spins_unit))
    return Structure(s, unit_cell, supercell=RotationSupercell(perpendicular, propagation_vector))

#
# Supercell methods
#

@check_sizes(directions=("n", 3), phases=("n",), force_numpy=True, allow_nones=True)
def propagation_vectors(
        directions: np.ndarray,
        phases: np.ndarray | None = None,
        incommensurate: bool=False) -> list[PropagationVector]:
    """Create list of propagation vectors.

    Parameters
    ----------
    directions : (n, 3) numpy.ndarray
        ``n``-by-3 matrix of directions (in lattice units).
    phases : (n,) numpy.ndarray, optional
        Length-``n`` vector of phases. Default is ``None``.
    incommensurate : bool, optional
        Whether or not to set them to be incommensurate. Default is ``False``.

    Returns
    -------
    list of PropagationVector
        Created propagation vectors.
    """
    if phases is None:
        phases = [0.0 for _ in range(directions.shape[0])]


    if incommensurate:
        constructor = PropagationVector
    else:
        constructor = CommensuratePropagationVector

    out = []
    for i in range(directions.shape[0]):
        pv = constructor(
                float(directions[i, 0]),
                float(directions[i, 1]),
                float(directions[i, 2]),
                phase=float(phases[i]))

        out.append(pv)

    return out

@check_sizes(perpendicular=(3,), propagation=(3,))
def helical_supercell(perpendicular: ArrayLike, propagation: ArrayLike):
    """Generate a helical supercell.

    Parameters
    ----------
    perpendicular : (3,) array_like
        Direction perpendicular to propagation vector, needed to specify "start"
        cell.
    propagation : (3,) array_like
        Propagation vector.

    Returns
    -------
    RotationSupercell
        Generated helical supercell.
    """
    return RotationSupercell(perpendicular=perpendicular, propagation_vector=propagation)

@check_sizes(directions=("n", 3), phases=("n",), axes=("n", 3), force_numpy=True, allow_nones=True)
def rotation_supercell(
        directions: ArrayLike,
        axes: ArrayLike,
        phases: ArrayLike | None = None,
        scaling: tuple[int, int, int]=(1,1,1)):
    """Create a supercell that rotates spins as we move along the propagation vector.

    Parameters
    ----------
    directions : (n, 3) array_like
        Propagation vector directions.
    axes : (n, 3) array_like
        Rotation axes for each propagation vector, or if only a single axis is
        specified, apply it to all of them.
    phases : (n,) array_like, optional
        Phases of the propagation vectors; ``0.0`` means starting with the spin as
        specified on the site. Default is ``None``.
    scaling : tuple of (int, int, int), optional
        Make a larger supercell by tiling the result this many times in each axis.
        Default is ``(1, 1, 1)``.

    Returns
    -------
    TransformationSupercell
        Supercell applying the specified spin rotations.
    """
    # Type of vectors will be assured to be list[CommensuratePropagationVector] as long as incommensurate is False
    vectors: list[CommensuratePropagationVector] = propagation_vectors(directions, phases, incommensurate=False)

    pairs = [(vector, RotationTransform(axes[i, :])) for i, vector in enumerate(vectors)]

    return TransformationSupercell(
        transforms=pairs,
        scaling=scaling)


@check_sizes(directions=("n", 3), phases=("n",), force_numpy=True, allow_nones=True)
def summation_supercell(
        directions: np.ndarray,
        phases: np.ndarray | None = None,
        scaling: tuple[int, int, int]=(1,1,1)):
    r"""Create a supercell based on the propagation vectors and partial spins.

    Specifically, :math:`m = \sum_j \mu_j \exp(2 \pi i d_j \cdot r + \phi_j)`.

    Parameters
    ----------
    directions : (n, 3) numpy.ndarray
        Propagation vector directions.
    phases : (n,) numpy.ndarray, optional
        Phases of the propagation vectors; ``0.0`` means starting with the spin as
        specified on the site. Default is ``None``.
    scaling : tuple of (int, int, int), optional
        Make a larger supercell by tiling the result this many times in each axis.
        Default is ``(1, 1, 1)``.

    Returns
    -------
    SummationSupercell
        Supercell summing partial spin components over the propagation vectors.
    """
    # Type of vectors will be assured to be list[CommensuratePropagationVector] as long as incommensurate is False
    vectors: list[CommensuratePropagationVector] = propagation_vectors(directions, phases, incommensurate=False)

    return SummationSupercell(vectors, scaling=scaling)

#
# Spacegroup methods
#

def get_spacegroup(search_string: str):
    """Get a spacegroup by name or symmetry operations.

    The searches are whitespace and case insensitive.

    Parameters
    ----------
    search_string : str
        Spacegroup symbol or semicolon-separated symmetry operations string.

    Returns
    -------
    SpaceGroup
        Matching spacegroup from the symmetry database.

    Raises
    ------
    NoSuchGroup
        If no matching spacegroup is found or multiple ambiguous matches exist.
    ValueError
        If the database lookup returns an unexpected match result type.

    Examples
    --------
    The following are equivalent::

        get_spacegroup("P1")
        get_spacegroup("p1")
        get_spacegroup("x,y,z")

    The following are equivalent::

        get_spacegroup("P-1")
        get_spacegroup("p-1")
        get_spacegroup("x,y,z; -x,-y,-z")

    Spacegroups with multiple settings sometimes need to have it specified::

        get_spacegroup("R3H")  # hexagonal setting
        get_spacegroup("R3R")  # rhombohedral setting
        get_spacegroup("B2/m")  # this setting of C2/m
        get_spacegroup("B 1 2/m 1")  # another way for the same group
        get_spacegroup("P 4/n 2/b 2/m : 1")  # setting 1
        get_spacegroup("P 4/n 2/b 2/m : 2")  # setting 2
    """
    if 'x' in search_string and 'y' in search_string and 'z' in search_string:
        # Not a spacegroup name, but might be a list of symmetry operations
        m = database.spacegroups_with_operations(search_string)

        if isinstance(m, ExactMatch):
            return m.spacegroup

        elif isinstance(m, PartialMatch):
            groups = m.spacegroups

            if len(groups) == 0:
                raise NoSuchGroup("Tried to find group by operations, no matching group found")

            elif len(groups) == 1:
                return groups[0]

            else:
                suggestions = ", ".join([group.symbol for group in groups[:-1]]) + " and " + groups[-1].symbol
                raise NoSuchGroup(f"Found {len(groups)} groups with these operations: {suggestions}")

        else:
            raise ValueError("Expected `spacegroups_by_operations` to return ExactMatch or PartialMatch")

    else:
        # Try to get by name
        return database.spacegroup_by_name(search_string)

@check_sizes(direction=(3,), force_numpy=True)
def filter(direction: ArrayLike,
           perpendicular: bool=False,
           symmetric: bool=False,
           max_dev_angle_deg: float=0.01) -> DirectionalityFilter:
    """Create a filter for directions (helper method for exchanges).

    Parameters
    ----------
    direction : (3,) array_like
        If not perpendicular, allowed direction of exchange. If perpendicular,
        normal to plane containing exchange.
    perpendicular : bool, optional
        Constrain to a line (``perpendicular=False``) or a plane
        (``perpendicular=True``). Default is ``False``.
    symmetric : bool, optional
        In the non-perpendicular case, ``symmetric=True`` generates exchanges in
        both directions. Default is ``False``.
    max_dev_angle_deg : float, optional
        Angular tolerance for the direction/normal in degrees. Default is ``0.01``.

    Returns
    -------
    DirectionalityFilter
        A ``DirectionalityFilter`` object that can be used to select exchanges in
        particular directions.
    """
    if perpendicular:
        return InPlaneFilter(direction=direction, max_dev_angle_deg=max_dev_angle_deg)
    else:
        if symmetric:
            return BiDirectionFilter(direction=direction, max_dev_angle_deg=max_dev_angle_deg)
        else:
            return InDirectionFilter(direction=direction, max_dev_angle_deg=max_dev_angle_deg)


def generate_exchanges(sites: list[Site] | Structure,
                       unit_cell: UnitCell | None = None,
                       exchange_type: type[Exchange] = HeisenbergExchange,
                       bond: int = 0,
                       max_distance: float = 0.0,
                       min_distance: float = 0.0,
                       direction_filter: DirectionalityFilter | None = None,
                       max_order: int | None = None,
                       j: float | None = None,
                       j_x: float | None = None,
                       j_y: float | None = None,
                       j_xy: float | None = None,
                       j_z: float | None = None,
                       d_x: float | None = None,
                       d_y: float | None = None,
                       d_z: float | None = None,
                       exchange_parameters: dict | None = None,
                       naming_pattern: str | None = None,
                       color: tuple[float, float, float] | None = None):
    """Automatically create a list of exchanges.

    Parameters
    ----------
    sites : list of Site or Structure
        List of sites to make exchanges between or a ``Structure`` object.
    unit_cell : UnitCell, optional
        Unit cell (needed if first argument is a list of sites). Default is ``None``.
    exchange_type : type of Exchange, optional
        Type of exchange. Default is ``HeisenbergExchange``.
    bond : int, optional
        The bond index (if this is given the ``*_distance`` parameters are ignored).
        Default is ``0``.
    max_distance : float, optional
        Maximum Cartesian distance (in Angstrom) at which exchanges are made.
        Default is ``0.0``.
    min_distance : float, optional
        Minimum Cartesian distance (in Angstrom) at which exchanges are made.
        Default is ``0.0``.
    direction_filter : DirectionalityFilter, optional
        Supply a ``DirectionalityFilter`` object (e.g. using ``filter``) to only
        create exchanges in certain directions. Default is ``None``.
    max_order : int, optional
        Maximum "order" of exchanges. Default is ``None``.
    j : float, optional
        Constant for scalar valued Heisenberg exchanges (only needed for
        Heisenberg). Default is ``None``.
    j_x : float, optional
        Constant for x component of Heisenberg-like exchanges (only needed for
        Diagonal). Default is ``None``.
    j_y : float, optional
        Constant for y component of Heisenberg-like exchanges (only needed for
        Diagonal). Default is ``None``.
    j_xy : float, optional
        Constant for x and y components of Heisenberg-like exchanges (only needed
        for XY and XXZ). Default is ``None``.
    j_z : float, optional
        Constant for z component of Heisenberg-like exchanges (only needed for
        Diagonal, Ising and XY). Default is ``None``.
    d_x : float, optional
        Constant for x component of Dzyaloshinskii-Moriya exchange (DM). Default is
        ``None``.
    d_y : float, optional
        Constant for y component of Dzyaloshinskii-Moriya exchange (DM). Default is
        ``None``.
    d_z : float, optional
        Constant for z component of Dzyaloshinskii-Moriya exchange (DM). Default is
        ``None``.
    exchange_parameters : dict, optional
        Parameters can be supplied in a dictionary instead. Default is ``None``.
    naming_pattern : str, optional
        String used to assign names to exchanges; see ``apply_naming_convention``.
        Default is ``None``.
    color : tuple of (float, float, float), optional
        Colour to draw the exchange with. Default is ``None``.

    Returns
    -------
    list of Exchange
        Generated list of exchanges.

    Raises
    ------
    RuntimeError
        If ``sites`` is not a ``Structure`` and ``unit_cell`` is ``None``, or if
        neither ``bond`` nor distance bounds are specified.
    """
    if isinstance(sites, Structure):
        unit_cell = sites.unit_cell
        sites = sites.sites
    elif unit_cell is None:
        raise RuntimeError('If first argument is not a Structure, you need to specify a unit_cell')

    if bond == 0 and max_distance == 0.0 and min_distance == 0.0:
        raise RuntimeError('You must specify either a bond index or maximum/minimum distance pairs')

    if exchange_parameters is None:
        exchange_parameters = {}

    if j is not None:
        exchange_parameters["j"] = j

    if j_x is not None:
        exchange_parameters["j_x"] = j_x

    if j_y is not None:
        exchange_parameters["j_y"] = j_y

    if j_z is not None:
        exchange_parameters["j_z"] = j_z

    if j_xy is not None:
        exchange_parameters["j_xy"] = j_xy

    if d_x is not None:
        exchange_parameters["d_x"] = d_x

    if d_y is not None:
        exchange_parameters["d_y"] = d_y

    if d_z is not None:
        exchange_parameters["d_z"] = d_z

    # Only want parameters that apply to the specific type
    used_parameters = {}
    for parameter, default in zip(exchange_type.parameters, exchange_type.parameter_defaults):
        if parameter in exchange_parameters:
            used_parameters[parameter] = exchange_parameters[parameter]
        else:
            used_parameters[parameter] = default

    # Adds extra parameters
    if color is not None:
        used_parameters['color'] = color

    group = ExchangeGroup(
        name = "<unnamed group>",
        bond = bond,
        min_distance = min_distance,
        max_distance = max_distance,
        max_order = max_order,
        naming_pattern = default_naming_pattern if naming_pattern is None else naming_pattern,
        exchange_type= exchange_type,
        exchange_parameters= used_parameters,
        direction_filter = direction_filter)

    return group.exchanges(sites, unit_cell)


@check_sizes(axis=(3, ), force_numpy=True)
def axis_anisotropies(
        sites: list[Site] | Structure,
        a: float,
        axis: ArrayLike = [0, 0, 1],
        name: str = ""):
    """Create anisotropy objects with magnitude ``a`` in direction ``axis`` for each site.

    Parameters
    ----------
    sites : list of Site or Structure
        Sites or ``Structure`` to attach anisotropies to.
    a : float
        Anisotropy magnitude.
    axis : (3,) array_like, optional
        Anisotropy axis direction. Default is ``[0, 0, 1]``.
    name : str, optional
        Name for the created anisotropy terms. Default is ``""``.

    Returns
    -------
    list of AxisMagnitudeAnisotropy
        Created axis-magnitude anisotropy objects for each site.
    """
    if isinstance(sites, Structure):
        sites = sites.sites
    return [AxisMagnitudeAnisotropy(site, a, axis, name) for site in sites]


@check_sizes(matrix=(3,3), force_numpy=True)
def matrix_anisotropies(
        sites: list[Site] | Structure,
        matrix: ArrayLike,
        name: str = ""):
    """Create anisotropy objects specified by a matrix, the same for each site.

    Parameters
    ----------
    sites : list of Site or Structure
        Sites or ``Structure`` to attach anisotropies to.
    matrix : (3, 3) array_like
        3x3 anisotropy matrix applied to each site.
    name : str, optional
        Name for the created anisotropy terms. Default is ``""``.

    Returns
    -------
    list of Anisotropy
        Created anisotropy objects for each site.
    """
    if isinstance(sites, Structure):
        sites = sites.sites
    return [Anisotropy(site, matrix, name) for site in sites]

def view(hamiltonian: Hamiltonian):
    """Show the current Hamiltonian in the viewer.

    Parameters
    ----------
    hamiltonian : Hamiltonian
        Hamiltonian object to display in the 3D viewer.
    """
    show_object(hamiltonian)

def load_tb2j(inputfile):
    """Load a TB2J output file and convert to a SpinW Hamiltonian object.

    Parameters
    ----------
    inputfile : str or os.PathLike
        Path to the TB2J output file.

    Returns
    -------
    Hamiltonian
        Converted ``Hamiltonian`` instance.
    """
    return TB2J_Input(inputfile).to_hamiltonian()
