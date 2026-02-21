""" Splitting and wrapping line segments to the unit cube"""

import numpy as np


def split_and_wrap_line_segment(a: np.ndarray, b: np.ndarray, tol=1e-12, max_parts=100_000):
    """ Split a line when it crosses integer values on the axes, and wrap it in 3D to [0,1]^3

    :return: (points where it crosses integers, line segments moved to [0,1])
    """
    a = np.array(a, dtype=float)
    b = np.array(b, dtype=float)

    current_point = a % 1 # Start in the box, regardless
    delta = b - a

    # We will need to be careful about zero change in a given direction,
    # if this is the case, we can ignore it and set it again later

    kept_dimensions = ~np.isclose(delta, 0)

    # Check for same point, need at least one dimension to work with
    if np.sum(kept_dimensions) == 0:
        return a.reshape(1, 3), [] # No line segments

    current_point = current_point[kept_dimensions]
    delta = delta[kept_dimensions]

    # We're going from one point to the other, so we'll always be going out one side, and in the other
    # but which one it is depends on the direction of delta
    # If the direction is positive, the next point in that dimension will be 1 and zero otherwise
    next_point_in_direction = 0.5*(1+np.sign(delta))

    pairs = []
    transition_points = [a[kept_dimensions]]

    # keep track of how far we need to go
    total_distance = np.abs(delta[0])

    for i in range(max_parts):
        # Step 1: Find distance in each component until the next integer

        dimension_distance = next_point_in_direction - current_point

        path_distance = dimension_distance / delta

        # Step 2: Find the smallest one
        #
        # masked = path_distance.copy()
        # masked[path_distance < tol] = np.inf
        dim = np.argmin(path_distance)

        # Step 3: Find the next point by moving path_distance * delta

        change = path_distance[dim] * delta
        next_transition_point = transition_points[-1] + change
        next_point = current_point + change

        # Step 4: Check for reaching the end of the line, two possibilities
        #  a) We have hit the end exactly (i.e. the line ends on an integer)
        #  b) We have overshot
        # Only check (b) now, as if we've hit it exactly, we want to save the next point

        if np.abs(next_transition_point[0] - transition_points[0][0]) >= total_distance + tol:
            b_kept_part = b[kept_dimensions]
            transition_points.append(b_kept_part)
            pairs.append((current_point, b_kept_part % 1))
            break


        # Step 5: Save, as long as the distance was not zero:
        if not np.any(path_distance < tol):
            transition_points.append(next_transition_point)
            pair = (current_point, next_point)
            pairs.append(pair)

        current_point = next_point.copy()

        # Step 6: Other exit state check (a)
        if np.isclose(np.abs(next_transition_point[0] - transition_points[0][0]), total_distance, atol=tol):
            # We've hit it dead on, and we've already saved everything, so just break
            break


        # Step 7: Shift, multiple dimensions can hit their target at the same time
        reached_target = np.isclose(current_point, next_point_in_direction, atol=tol)
        current_point[reached_target] = 1 - next_point_in_direction[reached_target]

    else:
        raise Exception("Iteration limit reached")

    # Finally, we need to re-introduce any dimensions we've ignored, can use either a or b for this, because
    #  we only keep dimensions where a=b
    transitions_array = np.tile(a.reshape(1, -1), (len(transition_points), 1))
    transitions_array[:, kept_dimensions] = np.array(transition_points)


    output_pairs = []
    for p, q in pairs:

        r = a.copy()
        s = b.copy()

        r[kept_dimensions] = q
        s[kept_dimensions] = p

        output_pairs.append((r, s))

    return transitions_array, output_pairs

if __name__ == "__main__":

    def run_example(a, b):
        """ Show plots and data for an example"""
        split, wrap = split_and_wrap_line_segment(a, b)

        import matplotlib.pyplot as plt

        print("Split points")
        print(split)
        print("Wrapped line segments")
        for pair in wrap:
            print(*pair)

        plt.subplot(1,2,1)
        plt.title("Split")
        for i in range(split.shape[0] - 1):
            plt.plot(split[i:i+2, 0], split[i:i+2, 1])

        plt.subplot(1, 2, 2)
        plt.title("Wrapped")
        for p, q in wrap:
            plt.plot([p[0], q[0]], [p[1], q[1]])

        plt.plot([0,0,1,1,0], [0,1,1,0,0], color='k') # box

        plt.show()


    run_example([0.5, 0.5, 0.5], [1.5, 2.5, 3.5]) # A nice check as it has a double collision at one point
    run_example([1.5, 2.5, 3.5], [0.5, 0.5, 0.5]) # Check reverse of it too


    run_example([0.1,0.5,0], [5, 0, 0]) # Lots of int values, one index doesn't change
    run_example([5, 0, 0], [0.1,0.5,0]) # And reverse

    run_example([-1, 1, -1], [1, -1, 1]) # Some negatives


