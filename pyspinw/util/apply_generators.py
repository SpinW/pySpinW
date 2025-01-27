import numpy as np


def _remove_duplicates_and_order_points(points: np.ndarray, tol=1e-8) -> np.ndarray:
    """ Assumes there are no duplicate sites with different moments"""

    # Remove duplicates first, as order might be very sensitive to rounding
    tol_squared = tol*tol

    # Step one, sort by coordinates

    points = points[points[:,0].argsort(), :] # Fastest search

    for i in range(1, 3):
        points = points[points[:, i].argsort(kind='mergesort'), :]  # mergesort keeps order from other dimensions

    # Step 2: collate points at same locations

    this_location = [points[0, :]]
    unique_locations = [this_location]
    for i in range(1, points.shape[0]):
        last_point = points[i-1, :]
        this_point = points[i, :]

        if np.sum((this_point[:3] - last_point[:3])**2) < tol_squared:
            this_location.append(this_point)
        else:
            this_location = [this_point]
            unique_locations.append(this_location)

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
                full_list.append(0.5*(location[0] + location[1]))
            case _:
                for point in location[1:]:
                    diff = np.sum((point - location[0]) ** 2)

                    if diff > tol_squared:
                        print(points)
                        print(location)
                        raise ValueError(f"Expected points with more than two copies to be exactly the same")

                full_list.append(location[0])

    # Return as an array
    return np.array(full_list, dtype=float)




def apply_generators_with_moments(points: np.ndarray, generators: list[tuple[np.ndarray, np.ndarray, float]], tol=1e-8) -> np.ndarray:

    for i, (quadratic, linear, time_reversal) in enumerate(generators):
        # print(f"Generator: {i+1} of {len(generators)}: {points.shape[0]}")
        # print(points)

        transformed_positions = points[:, :3] @ quadratic + linear
        transformed_positions %= 1 # Move back to unit cell

        transformed_moments = time_reversal * points[:, 3:]

        print(quadratic, linear, time_reversal)

        joined = np.concatenate((points,
                    np.concatenate((transformed_positions, transformed_moments), axis=1)), axis=0)

        points = _remove_duplicates_and_order_points(joined, tol=tol)

        # print(points)

    return points

def apply_generators_until_stable(
        points: np.ndarray,
        generators: list[tuple[np.ndarray, np.ndarray, float]], tol=1e-8,
        max_iters: int=1000) -> np.ndarray:

    points = _remove_duplicates_and_order_points(points)


    for i in range(max_iters): # Big but finite number
        # print(f"Repeated application: {i}")
        new_points = apply_generators_with_moments(points, generators)

        if new_points.shape == points.shape:

            difference = np.sum((points - new_points)**2)

            if difference < new_points.shape[0]*new_points.shape[1]*tol*tol:
                return points

        points = new_points

    raise Exception(f"Failed to reach steady state after {max_iters} iterations")

