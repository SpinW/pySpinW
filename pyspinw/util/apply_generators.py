import functools

import numpy as np

from pyspinw.checks import check_sizes
from pyspinw.util.group_generators import Generator

def _noisy_site_sort_comparator(tol):
    """ We want to sort a list of points so that equivalent points are next to each other,
    numerically, they only need to be within a tolerance, the sort comparison is therefore
    fairly non-standard

    We want to be able to provide a two place function with a fixed tolerance, therefore, it
    has been Curryed
    """

    def comparator(a: np.ndarray, b: np.ndarray):


        if np.abs(a[0] - b[0]) > tol:
            return 1 if a[0] < b[0] else -1
        else:
            if np.abs(a[1] - b[1]) > tol:
                return 1 if a[1] < b[1] else -1
            else:
                return 1 if a[2] < b[2] else -1

    return comparator

def _remove_duplicates_and_order_points(points: list[np.ndarray], tol) -> np.ndarray:
    """ Assumes there are no duplicate sites with different moments"""

    # Remove duplicates first, as order might be very sensitive to rounding
    tol_squared = tol*tol

    # Step one, sort by coordinates

    points = sorted(points, key=functools.cmp_to_key(_noisy_site_sort_comparator(tol)))
    # Step 2: collate points at same locations


    # import matplotlib.pyplot as plt
    # plt.subplot(1,2,1)
    # plt.scatter(points[:, 0], points[:, 1])
    # plt.subplot(1,2,2)
    # plt.scatter(points[:, 0], points[:, 2])
    # plt.show()

    last_point = points[0]
    this_location = [last_point]
    unique_locations = [this_location]

    for this_point in points[1:]:

        if np.sum((this_point[:3] - last_point[:3])**2) < tol_squared:
            this_location.append(this_point)
        else:
            this_location = [this_point]
            unique_locations.append(this_location)

        last_point = this_point

    # Go through the collections of points at each location,
    # and merge them in a way that makes their momentum zero if they are opposite
    # or the same if they are the same.

    # optional TODO: check that they are in fact the same up to sign

    full_list = []
    for location in unique_locations:
        match len(location):
            case 0:
                raise ValueError("Expected at least one point at each location but there's none, this should never happen")
            case 1:
                full_list += location
            case 2:
                if np.sum(np.abs(location[0][3:] - location[1][3:])) < tol:
                    full_list.append(location[0])
                elif np.sum(np.abs(location[0][3:] + location[1][3:])) < tol:
                    zeroed = location[0].copy()
                    zeroed[3:] = np.zeros((3,))
                    full_list.append(zeroed)
                else:
                    raise Exception("This is a problem")
            case _:
                print(location)
                raise ValueError(f"Expected at most two copies of a given site (with potentially opposing momenta)")

    # Return as an array
    return full_list





def _apply_generators_with_moments(
        points: list[np.ndarray],
        generators: list[Generator], tol) -> np.ndarray:

    for i, generator in enumerate(generators):
        # print(f"Generator: {i+1} of {len(generators)}: {points.shape[0]}")
        # print(points)

        transformed_positions = [generator(point.reshape(1, 6)).reshape(6) for point in points]

        joined = points + transformed_positions

        # print("\n\n\nPoints")
        # for point in points:
        #     print(point)
        #
        # print("\n\n\nJoined")
        # for point in joined:
        #     print(point)

        points = _remove_duplicates_and_order_points(joined, tol=tol)


    return points

@check_sizes(points=(-1, 6))
def apply_generators_with_moments(points: np.ndarray, generators: list[Generator], tol=1e-8) -> np.ndarray:
    point_list = [points[i, :] for i in range(points.shape[0])]
    point_list = _apply_generators_with_moments(point_list, generators, tol)

    return np.array(point_list)

def _apply_generators_until_stable(
        points: list[np.ndarray],
        generators: list[Generator],
        tol: float,
        max_iters: int) -> np.ndarray:


    for i in range(max_iters): # Big but finite number
        # print(f"Repeated application: {i}")
        new_points = _apply_generators_with_moments(points, generators, tol)

        if len(new_points) == len(points):
                return points

        points = new_points

    raise Exception(f"Failed to reach steady state after {max_iters} iterations")


@check_sizes(points=(-1, 6))
def apply_generators_until_stable(points: np.ndarray, generators: list[Generator], tol: float=1e-8, max_iters: int=1000) -> np.ndarray:
    point_list = [points[i, :] for i in range(points.shape[0])]
    point_list = _apply_generators_until_stable(point_list, generators, tol, max_iters)

    return np.array(point_list)