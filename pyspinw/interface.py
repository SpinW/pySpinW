""" Helper functions for python interface """
import math

import numpy as np
from numpy._typing import ArrayLike

from pyspinw.batch_couplings import default_naming_pattern
from pyspinw.checks import check_sizes
from pyspinw.coupling import Coupling
from pyspinw.couplinggroup import DirectionalityFilter, InPlaneFilter, InDirectionFilter, CouplingGroup
from pyspinw.site import LatticeSite
from pyspinw.symmetry.group import database, NoSuchGroup, ExactMatch, PartialMatch
from pyspinw.symmetry.supercell import PropagationVector, CommensuratePropagationVector, RotationTransform, \
    TransformationSupercell, SummationSupercell
from pyspinw.symmetry.unitcell import UnitCell


def sites(positions: ArrayLike,
          moments: ArrayLike | None = None,
          names: list[str] | None=None,
          convert_to_cell_with: UnitCell | None = None) -> list[LatticeSite]:

    """ Create lattice site

    :param positions: positions of the sites
    :param moments: moments of the sites, if not specified, they will be set to zero
    :param convert_to_cell_with: If this is None, we assume the position is in lattice units, and moments
                                are in the unit cell moment coordinate system
                                (see `pyspinw.UnitCell.moment_fractional_to_cartesian and
                                `pyspinw.UnitCell.moment_cartesian_to_fractional`)
                                if instead it is a unit cell, we assume the coordinates are cartesian,
                                and the positions will be converted into lattice units
    """

    # check the positions are n-by-3, or 1d
    positions = np.array(positions)
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
    moments = np.array(moments)
    match len(moments.shape):
        case 1:
            if n_sites != 1:
                raise ValueError("Expected same number of moments as positions")

            if moments.shape[0] != 3:
                raise ValueError("Expected moment to have three components")

            moments = moments.reshape(1,3,1)

        case 2:
            if n_sites != moments.shape[1]:
                raise ValueError("Expected same number of moments as positions")

            if moments.shape[1] != 3:
                raise ValueError("Expected moments to have three components")

            moments = moments.reshape(-1, 3, 1)

        case 3:
            if n_sites != moments.shape[1]:
                raise ValueError("Expected same number of moments as postions")

            if moments.shape[1] != 3:
                raise ValueError("Expected moments to have three components")

        case _:
            raise ValueError("Expected moments to be n-by-3-m, n-by-3, or vector of length 3")

    # Convert cell positions
    if convert_to_cell_with is not None:
        positions = convert_to_cell_with.cartesian_to_fractional(positions)

        # convert moments
        converted_moments = np.zeros_like(moments)
        for i in range(moments.shape[2]):
            converted_moments[:,:,i] = convert_to_cell_with.moment_cartesian_to_fractional(moments[:,:,i])

        moments = converted_moments

    n_sites = positions.shape[0]

    out = []
    for i in range(n_sites):
        position = positions[i, :]

        out.append(LatticeSite(
            i = position[0],
            j = position[1],
            k = position[2],
            supercell_moments=moments[i, :, :],
            name = names[i]
        ))

    return out

#
# Supercell methods
#

@check_sizes(directions=("n", 3), phases=("n",), force_numpy=True, allow_nones=True)
def propagation_vectors(
        directions: np.ndarray,
        phases: np.ndarray | None = None,
        incommensurate: bool=False) -> list[PropagationVector]:

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

    i.e. $m = \sum_j mu_j exp(2 \pi i d_j.r + \phi_j)$

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
def filter(direction: ArrayLike, perpendicular: bool=False, max_dev_angle_deg: float=0.01) -> DirectionalityFilter:
    if perpendicular:
        return InPlaneFilter(direction=direction, max_dev_angle_deg=max_dev_angle_deg)
    else:
        return InDirectionFilter(direction=direction, max_dev_angle_deg=max_dev_angle_deg)


def couplings(sites: list[LatticeSite],
              unit_cell: UnitCell,
              coupling_type: type[Coupling],
              max_distance: float,
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
              coupling_parameters: dict | None = None,
              naming_pattern: str | None = None,):

    """ Automatically creates a list of couplings

    :param sites: *required* List of sites to make couplings between
    :param unit_cell: *required* Unit cell (needed for working out the Cartesian coordinates of sites)
    :param coupling_type: *required* Type of coupling
    :param max_distance: *required* Maximum Cartesian distance at which couplings are made
    :param min_distance: Minimum Cartesian distance at which couplings are made
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
    :param coupling_parameters: Parameters can be supplied in a dictionary instead
    :param naming_pattern: String used to assign names to couplings, see `apply_naming_convention`

    :returns: list of couplings
    """

    if coupling_parameters is None:
        coupling_parameters = {}

    if j is not None:
        coupling_parameters["j"] = j

    if j_x is not None:
        coupling_parameters["j_x"] = j_x

    if j_y is not None:
        coupling_parameters["j_y"] = j_y

    if j_z is not None:
        coupling_parameters["j_z"] = j_z

    if j_xy is not None:
        coupling_parameters["j_xy"] = j_xy

    if d_x is not None:
        coupling_parameters["d_x"] = d_x

    if d_y is not None:
        coupling_parameters["d_y"] = d_y

    if d_z is not None:
        coupling_parameters["d_z"] = d_z

    # Only want parameters that apply to the specific type
    used_parameters = {}
    for parameter, default in zip(coupling_type.parameters, coupling_type.parameter_defaults):
        if parameter in coupling_parameters:
            used_parameters[parameter] = coupling_parameters[parameter]
        else:
            used_parameters[parameter] = default

    group = CouplingGroup(
        name = "<unnamed group>",
        min_distance = min_distance,
        max_distance = max_distance,
        max_order = max_order,
        naming_pattern = default_naming_pattern if naming_pattern is None else naming_pattern,
        coupling_type = coupling_type,
        coupling_parameters = used_parameters,
        direction_filter = direction_filter)

    return group.couplings(sites, unit_cell)