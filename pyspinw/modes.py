import numpy as np
from numpy._typing import ArrayLike

from pyspinw.path import Path

import numpy as np
from scipy.sparse import csr_matrix, coo_matrix
from scipy.sparse.csgraph import min_weight_full_bipartite_matching

# # sparse score matrix
# scores = csr_matrix([
#     [8, 0, 5],
#     [6, 7, 0],
#     [0, 9, 8]
# ])
#
# # convert max-score to min-cost
# cost = scores.max() - scores
#
# row_ind, col_ind = min_weight_full_bipartite_matching(cost)

# print(list(zip(row_ind, col_ind)))

def mode_separation_score(sorted_energies):
    """ Score for how well spread out energies are"""

    deltas = np.diff(sorted_energies)

    # Any power strictly greater than 2 will do,
    # needs delta to be positive unless it's an even power though
    return np.sum(deltas**2)

def _predictor(x_known: ArrayLike, y_known: ArrayLike, x_prediction: float):
    """ Predicts where a curve will go next based on its derivatives """

    x_known = np.array(x_known)
    y_known = np.array(y_known)

    match len(x_known):

        case 0:
            raise ValueError("Need at least one point")

        case 1:
            # Just one point, best prediction is that it is the same
            return y_known[-1]

        case 2:
            # Assume last dy/dx is the same
            dy_dx = (y_known[-1] - y_known[-2]) / (x_known[-1] - x_known[-2])
            delta_y = (x_prediction - x_known[-1]) * dy_dx
            return y_known[-1] + delta_y

        case 3:


            x_vals = x_known[-3:]
            y_vals = y_known[-3:]

            return (
            y_vals[0] * (x_prediction - x_vals[1]) * (x_prediction - x_vals[2]) /
            ((x_vals[0] - x_vals[1]) * (x_vals[0] - x_vals[2]))
            + y_vals[1] * (x_prediction - x_vals[0]) * (x_prediction - x_vals[2]) /
            ((x_vals[1] - x_vals[0]) * (x_vals[1] - x_vals[2]))
            + y_vals[2] * (x_prediction - x_vals[0]) * (x_prediction - x_vals[1]) /
            ((x_vals[2] - x_vals[0]) * (x_vals[2] - x_vals[1])))

        case _:
            x_vals = x_known[-4:]
            y_vals = y_known[-4:]

            return (
                    y_vals[0] * (x_prediction - x_vals[1]) * (x_prediction - x_vals[2]) * (x_prediction - x_vals[3]) /
                    ((x_vals[0] - x_vals[1]) * (x_vals[0] - x_vals[2]) * (x_vals[0] - x_vals[3]))

                    + y_vals[1] * (x_prediction - x_vals[0]) * (x_prediction - x_vals[2]) * (x_prediction - x_vals[3]) /
                    ((x_vals[1] - x_vals[0]) * (x_vals[1] - x_vals[2]) * (x_vals[1] - x_vals[3]))

                    + y_vals[2] * (x_prediction - x_vals[0]) * (x_prediction - x_vals[1]) * (x_prediction - x_vals[3]) /
                    ((x_vals[2] - x_vals[0]) * (x_vals[2] - x_vals[1]) * (x_vals[2] - x_vals[3]))

                    + y_vals[3] * (x_prediction - x_vals[0]) * (x_prediction - x_vals[1]) * (x_prediction - x_vals[2]) /
                    ((x_vals[3] - x_vals[0]) * (x_vals[3] - x_vals[1]) * (x_vals[3] - x_vals[2]))
            )



def mode_sort(path: Path, energies: list[np.ndarray], intensities: list[np.ndarray], cutoff_ratio=2.0):
    """ Try to make continuous modes """

    # number of modes
    n_modes = len(energies[0])

    # Q values
    q_points = path.q_points()
    delta_q = q_points[1:, :] - q_points[:-1, :]
    q_distances = np.sqrt(np.sum(delta_q**2, axis=1))
    path_length = np.concatenate(([0.0], q_distances))

    # Output
    output_energies = [[] for _ in range(n_modes)]
    output_intensities = [[] for _ in range(n_modes)]


    # We can only expect modes to be continuous over a single direction in q space, so
    # we split into sections based on the direction of the path
    # It is safe to assume that the path points will be used to define unique directions
    for slice in path.section_slices():
        print(slice)

        q = path_length[slice]
        q -= q[0]

        section_energies = energies[slice]
        section_intensities = intensities[slice]

        energy_intensity_pairs = [np.array([energy, intensity])
                                  for energy, intensity in zip(section_energies, section_intensities)]

        # Step 1, for sparseness, we need to get an estimate of what a sensible cutoff is
        # To do this, make an assumption that the q resolution is high enough that for most
        # of the time, the sorted energies will be close to what we want except for at certain
        # crossing points, this way we can measure the gradient of the modes and get a sensible
        # maximum change per q division

        # (n_q, 2, n_modes)
        sorted_pairs = np.array([pair[:, np.argsort(pair[0, :])] for pair in energy_intensity_pairs])
        energy_runs = sorted_pairs[:, 0, :]
        intensity_runs = sorted_pairs[:, 1, :]

        diffs = np.abs(energy_runs[1:, :] - energy_runs[:-1, :])
        cutoff = cutoff_ratio*np.max(diffs)

        # We want to find a good starting point, one where they're as spread out in energy as possible

        scores = [mode_separation_score(pair[0, :]) for pair in sorted_pairs]
        starting_point = np.argmin(scores)

        # Degenerate modes need to have some factor that keeps them having the right intensity

        starting_energy = energy_runs[starting_point, :]
        starting_intensity = intensity_runs[starting_point, :]

        mode_qs = [q[starting_point]]
        mode_energies = [[energy] for energy in starting_energy]
        mode_intensities = [[intensity] for intensity in starting_intensity]

        for partial_range, go_forwards in [
            (range(starting_point + 1, len(q)), True),
            (range(starting_point - 1, -1, -1), False)]:

            print(partial_range, go_forwards)

            step = 1 if go_forwards else -1

            for i in partial_range:

                # Predict the next energy and intensity values for each mode
                predicted_energies = [_predictor(mode_qs, series[::step], q[i]) for series in mode_energies]
                predicted_intensities = [_predictor(mode_qs, series[::step], q[i]) for series in mode_intensities]

                # Score them by distance from actual values

                energy_differences = np.abs(np.array(predicted_energies)[None, :] - sorted_pairs[i, 0, :, None])
                intensity_differences = np.abs(np.array(predicted_intensities)[None, :] - sorted_pairs[i, 1, :, None])

                small_enough = energy_differences < cutoff
                match_matrix = np.zeros_like(energy_differences)

                match_matrix[small_enough] = energy_differences[small_enough] + 1
                # match_matrix[~small_enough] = float('inf')

                # Intensity differences will need to be accounted for by some correction

                # print(energy_differences)
                # print(match_matrix)

                match_matrix_coo = coo_matrix(match_matrix.T)

                #print(match_matrix_coo)

                row_ind, col_ind = min_weight_full_bipartite_matching(match_matrix_coo)

                #print(row_ind, col_ind)

                # Now, apply the sorting
                if go_forwards:
                    for row, col in zip(row_ind, col_ind):
                        mode_energies[row].append(sorted_pairs[i, 0, col])
                        mode_intensities[row].append(sorted_pairs[i, 1, col])

                else:
                    for row, col in zip(row_ind, col_ind):
                        mode_energies[row].insert(0, sorted_pairs[i, 0, col])
                        mode_intensities[row].insert(0, sorted_pairs[i, 1, col])

        # join section modes together

        for i in range(n_modes):
            output_energies[i] += mode_energies[i]
            output_intensities[i] += mode_intensities[i]

    return output_energies, output_intensities



if __name__ == "__main__":

    predictor_check = True
    assignment_check = True

    if predictor_check:

        import matplotlib.pyplot as plt

        n = 20
        x = np.linspace(0, 2*np.pi, n)
        y = np.sin(x)

        plt.plot(x, y)
        x_pred = []
        y_pred = []
        for i in range(1,n):
            x_pred.append(x[i])
            predicted_y = _predictor(x[:i], y[:i], x[i])
            y_pred.append(predicted_y)

        plt.scatter(x_pred, y_pred)
        plt.show()


    if assignment_check:
        path = Path([[0,0,0],[0,0,1],[0,1,1]])

        points = path.q_points()

        x = points[:,1]
        y = points[:,2]

        a = np.sin(2*np.pi*x) + np.sin(4*np.pi*y) + 1
        b = np.sin(4*np.pi*x) + 1

        energies = [[aa, bb, -bb, -aa] for aa, bb in zip(a, b)]
        intensities = [[1, 1, 1, 1] for aa, bb in zip(a, b)]

        energies, intensities = mode_sort(path, energies, intensities)


        import matplotlib.pyplot as plt
        plt.figure("Input data")

        plt.plot(path.x_values(), a)
        plt.plot(path.x_values(), b)
        plt.plot(path.x_values(), -b)
        plt.plot(path.x_values(), -a)

        plt.figure("Output data")

        for mode_energy in energies:
            plt.plot(mode_energy)

        plt.show()

