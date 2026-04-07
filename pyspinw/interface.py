""" Helper functions for python interface """
import math

import numpy as np
from numpy._typing import ArrayLike

from pyspinw.anisotropy import AxisMagnitudeAnisotropy, Anisotropy
from pyspinw.batch_couplings import default_naming_pattern
from pyspinw.checks import check_sizes
from pyspinw.exchange import Exchange, HeisenbergExchange
from pyspinw.exchangegroup import DirectionalityFilter, InPlaneFilter, InDirectionFilter, ExchangeGroup, \
    BiDirectionFilter
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.site import LatticeSite
from pyspinw.structures import Structure
from pyspinw.symmetry.group import database, NoSuchGroup, ExactMatch, PartialMatch
from pyspinw.symmetry.supercell import PropagationVector, CommensuratePropagationVector, RotationTransform, \
    TransformationSupercell, SummationSupercell, RotationSupercell
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.units import CoordsUnits

try:
    from pyspinw.gui.viewer import show_hamiltonian
except ModuleNotFoundError:
    def show_hamiltonian(*args, **kwargs):
        """ Show a Hamiltonian in the viewer"""
        raise RuntimeError('PySide or OpenGL not installed')

def _check_positions_moments_shape(positions: ArrayLike,
                                   moments: ArrayLike | None = None) -> tuple[ArrayLike, ArrayLike]:
    """ Checks position and moment arrays and put them into canonical shape """
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

    # check moment definitions
    if moments is None:
        moments = np.zeros((n_sites, 3, 1))
    else:
        moments = np.array(moments, dtype=float)
    match len(moments.shape):
        case 1:
            if n_sites != 1:
                raise ValueError("Expected same number of moments as positions")

            if moments.shape[0] != 3:
                raise ValueError("Expected moment to have three components")

            moments = moments.reshape(1,3,1)

        case 2:
            if n_sites != moments.shape[0]:
                raise ValueError("Expected same number of moments as positions")

            if moments.shape[1] != 3:
                raise ValueError("Expected moments to have three components")

            moments = moments.reshape(-1, 1, 3)

        case 3:
            if n_sites != moments.shape[0]:
                raise ValueError("Expected same number of moments as postions")

            if moments.shape[2] != 3:
                raise ValueError("Expected moments to have three components")

        case _:
            raise ValueError("Expected moments to be n-by-3-m, n-by-3, or vector of length 3")

    return positions, moments

def generate_sites(positions: ArrayLike,
                   moments: ArrayLike | None = None,
                   names: list[str] | None=None) -> list[LatticeSite]:
    """ Create a list of lattice sites

    :param positions: positions of the sites
    :param moments: moments of the sites, if not specified, they will be set to zero
    :param names: (optional) a list of names for sites
    """
    positions, moments = _check_positions_moments_shape(positions, moments)

    n_sites = positions.shape[0]
    if names is None:
        names = [f'site{i:02d}' for i in range(n_sites)]

    out = []
    for i in range(n_sites):
        position = positions[i, :]

        out.append(LatticeSite(
            i = position[0],
            j = position[1],
            k = position[2],
            supercell_moments=np.squeeze(moments[i, :, :]),
            name = names[i]
        ))

    return out

def _transform_site(unit_cell: UnitCell,
                    positions: ArrayLike,
                    moments: ArrayLike,
                    names: list[str] | None,
                    magnitudes: ArrayLike | None,
                    positions_unit: CoordsUnits | str,
                    moments_unit: CoordsUnits | str) -> tuple:
    """ Converts site data into default units """
    positions, moments = _check_positions_moments_shape(positions, moments)
    # check units
    positions_unit = CoordsUnits(positions_unit)
    moments_unit = CoordsUnits(moments_unit)
    # Convert cell positions
    if positions_unit == CoordsUnits.XYZ:
        positions = unit_cell.cartesian_to_fractional(np.array(positions))
    # convert moments
    if moments_unit == CoordsUnits.LU:
        moments = np.array(moments, dtype=float)
        for i in range(moments.shape[0]):
            magnitude = np.linalg.norm(moments[i,:,:]) if magnitudes is None else magnitudes[i]
            moment = unit_cell.fractional_to_cartesian(moments[i,:,:])
            moments[i,:,:] = moment * (magnitude / np.linalg.norm(moment))
    elif magnitudes is not None:
        for i in range(moments.shape[0]):
            moments[i,:,:] *= (magnitudes[i] / np.linalg.norm(moments[i,:,:]))
    return positions, moments, names

def generate_structure(unit_cell: UnitCell,
                       positions: ArrayLike,
                       moments: ArrayLike,
                       propagation_vectors: ArrayLike | None = None,
                       names: list[str] | None = None,
                       magnitudes: ArrayLike | None = None,
                       positions_unit: CoordsUnits | str = 'lu',
                       moments_unit: CoordsUnits | str = 'xyz') -> Structure:
    """ Creates a magnetic structure from a list of positions and moments

    :param unit_cell: a UnitCell object
    :param positions: positions of the sites
    :param moments: moments of the sites
    :param perpendicular: vector perpendicular to helix plane of rotation
    :param propagation_vectors: the propagation vector(s)
    :param names: (optional) a list of names for sites
    :param magnitudes: (optional) a list of moment magnitudes (spin length S)
                       If not specified, the norm of the moments vectors will be used as the spin length
    :param positions_units: the units for atomic positions; either 'lu' or 'xyz'
                            (default: 'lu', lattice units)
    :param moments_units: the units for the magnetic moment directions; either 'lu' or 'xyz'
                          (default: 'xyz', Cartesian axis with x||a)
    :param convert_to_cell_with: (optional) Must be specified if positions_units is not 'lu' and
                                 moments_units is not 'xyz' and allows conversion from those units
                                 to the default used in LatticeSite.
                                 (see `pyspinw.UnitCell.moment_fractional_to_cartesian and
                                 `pyspinw.UnitCell.moment_cartesian_to_fractional`)
    """
    s = generate_sites(*_transform_site(unit_cell, positions, moments, names, magnitudes, positions_unit, moments_unit))
    if propagation_vectors is not None:
        pvs = np.array(propagation_vectors, ndmin=2)
        if pvs.shape[0] != 3 and pvs.shape[1] == 3:
            pvs = pvs.T
        pvs = [CommensuratePropagationVector(*pvs[:,ip]) for ip in range(pvs.shape[1])]
        return Structure(s, unit_cell, supercell=SummationSupercell(propagation_vectors=pvs))
    return Structure(s, unit_cell)

def generate_helical_structure(unit_cell: UnitCell,
                               positions: ArrayLike,
                               moments: ArrayLike,
                               perpendicular: ArrayLike,
                               propagation_vector: ArrayLike,
                               names: list[str] | None=None,
                               magnitudes: ArrayLike | None = None,
                               positions_unit: CoordsUnits | str = 'lu',
                               moments_unit: CoordsUnits | str = 'xyz') -> Structure:
    """ Creates a helical structure with a propagation vector and plane normal

    :param unit_cell: a UnitCell object
    :param positions: positions of the sites
    :param moments: moments of the sites
    :param perpendicular: vector perpendicular to helix plane of rotation
    :param propagation_vector: the propagation vector
    :param names: (optional) a list of names for sites
    :param magnitudes: (optional) a list of moment magnitudes (spin length S)
                       If not specified, the norm of the moments vectors will be used as the spin length
    :param positions_units: the units for atomic positions; either 'lu' or 'xyz'
                            (default: 'lu', lattice units)
    :param moments_units: the units for the magnetic moment directions; either 'lu' or 'xyz'
                          (default: 'xyz', Cartesian axis with x||a)
    :param convert_to_cell_with: (optional) Must be specified if positions_units is not 'lu' and
                                 moments_units is not 'xyz' and allows conversion from those units
                                 to the default used in LatticeSite.
                                 (see `pyspinw.UnitCell.moment_fractional_to_cartesian and
                                 `pyspinw.UnitCell.moment_cartesian_to_fractional`)
    """
    s = generate_sites(*_transform_site(unit_cell, positions, moments, names, magnitudes, positions_unit, moments_unit))
    return Structure(s, unit_cell, supercell=RotationSupercell(perpendicular, propagation_vector))

#
# Supercell methods
#

@check_sizes(directions=("n", 3), phases=("n",), force_numpy=True, allow_nones=True)
def propagation_vectors(
        directions: np.ndarray,
        phases: np.ndarray | None = None,
        incommensurate: bool=False) -> list[PropagationVector]:
    """ Create list of propagation vectors

    :param directions: n-by-3 matrix of directions (in lattice units)
    :param phases: n vector of phases
    :param incommensurate: Whether or not to set them to be incommensuate
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

@check_sizes(perpendicular=(3,), propagation=(3,))
def helical_supercell(perpendicular: ArrayLike, propagation: ArrayLike):
    """ Generate a helical supercell

    :param perpendicular: (3-vector) Direction perpendicular to propagation vector, needed to specify "start" cell
    :param propagation: (3-vector) Propagation vector

    """
    return RotationSupercell(perpendicular=perpendicular, propagation_vector=propagation)

@check_sizes(directions=("n", 3), phases=("n",), force_numpy=True, allow_nones=True)
def rotation_supercell(
        directions: ArrayLike,
        axes: ArrayLike,
        phases: ArrayLike | None = None,
        scaling: tuple[int, int, int]=(1,1,1)):
    """ Create a supercell that rotates spins as we move along the propagation vector

    :param directions: propagation vector directions
    :param axes: rotation axes for each propagation vector,
                 or if only a single axis is specified, apply it to all of them
    :param phases: Phases of the propagation vectors, 0.0 means starting with the moment as specified on the site
    :param scaling: Make a larger supercell by tiling the result this many times in each axis
    """
    # Check that the axes match up

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
    """ Create a supercell based on the propagation vectors and partial moments

    i.e. $m = \\sum_j mu_j exp(2 \\pi i d_j.r + \\phi_j)$

    :param directions: propagation vector directions
    :param phases: Phases of the propagation vectors, 0.0 means starting with the moment as specified on the site
    :param scaling: Make a larger supercell by tiling the result this many times in each axis
    """
    # Type of vectors will be assured to be list[CommensuratePropagationVector] as long as incommensurate is False
    vectors: list[CommensuratePropagationVector] = propagation_vectors(directions, phases, incommensurate=False)

    return SummationSupercell(vectors, scaling=scaling)

#
# Spacegroup methods
#

def spacegroup(search_string: str):
    """ Get a spacegroup by name or symmetry operations

    The searches are whitespace and case insensitive.

    Examples:
        The following are equivalent:
            spacegroup("P1")
            spacegroup("p1")
            spacegroup("x,y,z")

        The following are equivalent:
            spacegroup("P-1")
            spacegroup("p-1")
            spacegroup("x,y,z; -x,-y,-z")

        Spacegroups with multiple settings sometimes need to have it specified, but some don't:
            spacegroup("R3")  **fails**
            spacegroup("R3H")  **hexagonal setting**
            spacegroup("R3R")  **rhombohedral setting**

            spacegroup("B2/m") **this setting of C2/m**
            spacegroup("B 1 2/m 1") **another way for the same group**

            spacegroup("P 4/n 2/b 2/m : 1") **setting 1**
            spacegroup("P 4/n 2/b 2/m : 1") **setting 2**

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

            elif 1 < len(groups) <= 5:
                suggestions = ", ".join([group.symbol for group in groups[:-1]]) + " and " + groups[-1].symbol
                raise NoSuchGroup(f"Found multiple groups with these operations: {suggestions}")

            else:
                raise NoSuchGroup(f"Found {len(groups)} groups matching those operations.")

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
    """ Create a filter for directions (helper method for couplings)

    :param direction: If not perpendicular, allowed direction of coupling
                      If perpendicular, normal to plane containing coupling
    :param perpendicular: Constrain to a line (perpendicular=False) or a plane (perpendicular=True)
    :param symmetric: In the not perpendicular case, symmetric True generates couplings in both directions
    :param max_dev_angle_deg: Angular tolerance for the direction/normal in degrees

    :returns: A DirectionalityFilter object that can be used to select couplings in particular directions
    """
    if perpendicular:
        return InPlaneFilter(direction=direction, max_dev_angle_deg=max_dev_angle_deg)
    else:
        if symmetric:
            return BiDirectionFilter(direction=direction, max_dev_angle_deg=max_dev_angle_deg)
        else:
            return InDirectionFilter(direction=direction, max_dev_angle_deg=max_dev_angle_deg)


def generate_exchanges(sites: list[LatticeSite] | Structure,
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
                       naming_pattern: str | None = None, ):
    """ Automatically creates a list of couplings

    :param sites: *required* List of sites to make couplings between or a Structure object
    :param unit_cell: Unit cell (needed if first argument is a list of sites)
    :param exchange_type: Type of coupling (defaults to HeisenbergCoupling)
    :param bond: The bond index (If this is given the _distance parameters are ignored)
    :param max_distance: Maximum Cartesian distance (in Angstrom) at which couplings are made
    :param min_distance: Minimum Cartesian distance (in Angstrom) at which couplings are made
    :param direction_filter: Supply a DirectionalityFilter object (e.g. using `filter`)
                             to only create couplings in certain directions
    :param max_order: Maximum "order" of couplings
    :param j" Constant for scalar valued Heisenberg couplings (only needed for Heisenberg)
    :param j_x" Constant for x component of Heisenberg-like couplings (only needed for Diagonal)
    :param j_y" Constant for y component of Heisenberg-like couplings (only needed for Diagonal)
    :param j_z" Constant for z component of Heisenberg-like couplings (only needed for Diagonal, Ising and XY)
    :param j_xy" Constant for x and y components of Heisenberg-like couplings (only needed for XY and XXZ)
    :param d_x" Constant for x component of Dzyaloshinskii-Moriya coupling (DM)
    :param d_y" Constant for y component of Dzyaloshinskii-Moriya coupling (DM)
    :param d_z" Constant for z component of Dzyaloshinskii-Moriya coupling (DM)
    :param exchange_parameters: Parameters can be supplied in a dictionary instead
    :param naming_pattern: String used to assign names to couplings, see `apply_naming_convention`

    :returns: list of couplings
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

    group = ExchangeGroup(
        name = "<unnamed group>",
        bond = bond,
        min_distance = min_distance,
        max_distance = max_distance,
        max_order = max_order,
        naming_pattern = default_naming_pattern if naming_pattern is None else naming_pattern,
        exchange_type= exchange_type,
        coupling_parameters = used_parameters,
        direction_filter = direction_filter)

    return group.exchanges(sites, unit_cell)


@check_sizes(axis=(3, ), force_numpy=True)
def axis_anisotropies(
        sites: list[LatticeSite] | Structure,
        a: float,
        axis: ArrayLike = [0, 0, 1]):
    """ Create anisotropy objects with magnitude `a` in direction `axis` for each site """
    if isinstance(sites, Structure):
        sites = sites.sites
    return [AxisMagnitudeAnisotropy(site, a, axis) for site in sites]


@check_sizes(matrix=(3,3), force_numpy=True)
def matrix_anisotropies(
        sites: list[LatticeSite] | Structure,
        matrix: ArrayLike):
    """ Create anisotropy objects specified by a matrix, the same for each site """
    if isinstance(sites, Structure):
        sites = sites.sites
    return [Anisotropy(site, matrix) for site in sites]

def view(hamiltonian: Hamiltonian):
    """ Show the current Hamiltonian in the viewer"""
    show_hamiltonian(hamiltonian)
