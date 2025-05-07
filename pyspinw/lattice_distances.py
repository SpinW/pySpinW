"""Calculations for finding distances between sites in a lattice"""

from dataclasses import dataclass

import numpy as np
from pyspinw.checks import check_sizes

@dataclass
class InteractionGeometries:
    """Holds the output from the algorithm to find relative positions"""

    cell_indices: np.ndarray
    """ Indices of the cells where a point was found at a sufficiently close distance"""

    vectors: np.ndarray
    """ Vectors for each position found"""

    distances: np.ndarray
    """ Distance to the each point"""


@check_sizes(fractional_coordinates=(3,),
             unit_cell_transform=(3,3))
def find_relative_positions(
        fractional_coordinates: np.ndarray,
        unit_cell_transform: np.ndarray,
        max_distance: float,
        tol: float=1e-7,
        allow_self=True) -> InteractionGeometries:
    """Find fractional coordinates from translations of the unit cell close to the origin.

    Find all sets of fractional coordinates produced by translations of input coordinates by multiples
    of the unit cell, and that differ by a cartesian distance of a most max_distance from (0,0,0)

    :param fractional_coordinates: Fractional coordinates of the base point
    :param unit_cell_transform: Transform from fractional coordinates to cartesian coordinates
    :param max_distance: Maximum distance
    :param tol: tolerance used to check for identity
    :param allow_self: include the original point with no translation (i,j,k = 0,0,0)
    """
    fractional_coordinate_offsets = get_cell_offsets_containing_bounding_box(unit_cell_transform, max_distance)
    fractional_positions = fractional_coordinate_offsets + fractional_coordinates.reshape(1, 3)

    cartesian_position = fractional_positions @ unit_cell_transform.T

    square_distances = np.sum(cartesian_position**2, axis=1)

    within_distance = square_distances <= (max_distance + tol)**2

    if allow_self:
        valid_distance = within_distance
    else:
        not_self = square_distances > tol*tol
        valid_distance = np.logical_and(within_distance, not_self)


    output_indices = np.array(fractional_coordinate_offsets[valid_distance, :], dtype=int)
    output_positions = cartesian_position[valid_distance, :]
    output_distances = np.sqrt(square_distances[valid_distance])

    return InteractionGeometries(output_indices, output_positions, output_distances)

def get_cell_offsets_containing_bounding_box(
        unit_cell_transform: np.ndarray,
        radius: np.ndarray) -> np.ndarray:
    """List of unit cell translations that corresponds to

    :param unit_cell_transform: transformation from fractional coordinates to cartesian coordinates
    :param radius: maximum radius
    """
    # Note: In fractional coordinates, a fixed radius sphere is an ellipsoid
    # Note: We want any unit cell that overlaps with the ellipsoid, as such, we need to add
    #       a value of 0.5 at many points

    # The ellipse is Ax^2 + By^2 + Cz^2 + 2Dxy + 2Eyz + 2Fxz = 0
    # i.e. if T is the transform such that
    #   T.(i,j,k) = (x,y,z)
    # and V is the vector (x,y,z) then we have
    #   (TV).(TV) = 0
    # or equivalently
    #   V.MV = 0   with M = T.T
    # and M is
    #   [[A D F]
    #    [D B E]
    #    [F E C]]
    #

    # Get the maximum values to search in each cartesian direction
    # The bounding box for the ellipsoid is given by
    #  V_max[i] = r sqrt((Q^-1)[i])
    # we need to add 1/2 to this before we do a search
    #

    box_sizes = (radius *
                 np.sqrt(
                     np.diag(
                        np.linalg.inv(
                            np.dot(
                                unit_cell_transform,
                                unit_cell_transform)))))

    # Build a grid of points within the bounding box: we want the corners of the boxes
    # The first box is centred on 0,0,0, and we want at least one on point for each box outside
    #  the ellipse at the extremes
    #

    limits = np.ceil(box_sizes)

    # Get bounds of the form [-(l+1), l+1], note that arange is not inclusive
    i_values = np.arange(-(limits[0]+1), limits[0]+2)
    j_values = np.arange(-(limits[1]+1), limits[1]+2)
    k_values = np.arange(-(limits[2]+1), limits[2]+2)

    i, j, k = np.meshgrid(i_values, j_values, k_values)

    ijk = np.concatenate((
        i.reshape(-1,1),
        j.reshape(-1,1),
        k.reshape(-1,1)
    ), axis=1)

    return ijk


def demo_point_finding():
    """Demonstrate the calculation of the calculation of cells for searching in"""
    # We want to only import matplotlib here
    #pylint: disable=C0415
    import matplotlib.pyplot as plt
    # pylint: enable=C0415

    radius = 5
    # transform =np.array(object=
    #     [[ 2, 0,          0         ],
    #      [ 1, 1.73205078, 0         ],
    #      [ 1, 0.57735026, 1.63299322]])

    transform =np.array(object=
        [[ 1, -0.1, 0.02 ],
         [ 0.6, 3, 0.02 ],
         [ 1, 1.3, 1.2 ]])

    point_fractional_position = np.array([0.2, 0.4, 0.6])
    # point_fractional_position = np.array([0, 0, 0], dtype=float)

    #
    # Plot a slice of the ellipse at z=0 in fractional coordinates
    #

    angles = np.linspace(0, 2*np.pi, 101).reshape(1, -1)
    untransformed_slice = radius * np.concatenate(
        (np.sin(angles),
         np.cos(angles),
         np.zeros_like(angles)), axis=0)

    transformed_limits = np.linalg.inv(transform) @ untransformed_slice

    plt.figure("Lattice cells checked")
    plt.plot(transformed_limits[0, :], transformed_limits[1, :])



    #
    # Plot the initial cell corner points used for searching
    #

    initial_points = get_cell_offsets_containing_bounding_box(transform, radius)

    # Take a slice for drawing
    z_plane = initial_points[np.floor(initial_points[:, 2]) == 0, :]
    plt.scatter(z_plane[:, 0], z_plane[:, 1])

    # Find the relative positions
    #

    results = find_relative_positions(
        fractional_coordinates=point_fractional_position,
        unit_cell_transform=transform,
        max_distance=radius)

    # Plot the indices that result in sufficiently short distance
    inds = results.cell_indices
    plt.scatter(inds[:, 0], inds[:, 1])

    #

    plt.figure("Cartesian Space")

    # Plot the cicle
    plt.plot(untransformed_slice[0, :], untransformed_slice[1, :])

    # Plot the checked points
    fractional_positions = initial_points + point_fractional_position.reshape(1, 3)
    cartesian_position = fractional_positions @ transform.T

    plt.scatter(cartesian_position[:, 0], cartesian_position[:, 1])

    # Plot the chosen points
    plt.scatter(results.vectors[:, 0], results.vectors[:, 1])

    plt.show()

if __name__ == "__main__":
    demo_point_finding()
