import numpy as np


def _remove_duplicates_and_order_points(points: np.ndarray, tol=1e-8) -> np.ndarray:
    # Remove duplicates first, as order might be very sensitive to rounding
    tol_squared = tol*tol

    unique_points = []
    for i in range(points.shape[0]):
        test_point = points[i, :]
        for unique_point in unique_points:
            if np.sum((unique_point - test_point)**2) < tol_squared:
                break
        else:
            unique_points.append(test_point)

    output = np.array(unique_points)
    output = output[output[:,0].argsort(), :] # Fastest search

    for i in range(1, output.shape[1]):
        output = output[output[:, i].argsort(kind='mergesort'), :]  # mergesort keeps order from other dimensions

    return output



def apply_generators_with_moments(points: np.ndarray, generators: list[tuple[np.ndarray, np.ndarray, float]]) -> np.ndarray:

    for quadratic, linear, time_reversal in generators:

        transformed_positions = points[:, :3] @ quadratic + linear
        transformed_positions %= 1 # Move back to unit cell

        transformed_moments = time_reversal * points[:, 3:]


        joined = np.concatenate((points,
                    np.concatenate((transformed_positions, transformed_moments), axis=1)), axis=0)

        points = _remove_duplicates_and_order_points(joined)

    return points

def apply_generators_until_stable(
        points: np.ndarray,
        generators: list[tuple[np.ndarray, np.ndarray, float]], tol=1e-8,
        max_iters: int=1000) -> np.ndarray:

    points = _remove_duplicates_and_order_points(points)


    for i in range(max_iters): # Big but finite number
        # print(i)
        new_points = apply_generators_with_moments(points, generators)

        if new_points.shape != points.shape:
            points = new_points
            continue

        difference = np.sum((points - new_points)**2)

        if difference < new_points.shape[0]*new_points.shape[1]*tol*tol:
            return points

    raise Exception(f"Failed to reach steady state after {max_iters} iterations")

