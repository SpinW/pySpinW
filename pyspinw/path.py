import numpy as np
from numpy._typing import ArrayLike

from pyspinw.checks import check_sizes


class Path:
    """ Path through q-space"""

    def __init__(self,
                 points: ArrayLike,
                 labels: list[str] | None = None,
                 avoid_endpoints=True,
                 resolution: int=101,
                 scale_by_distance=False):

        self._points = np.array(points, dtype=float)

        if len(self._points.shape) != 2 or self._points.shape[1] != 3:
            raise ValueError("Expected points to be a list of 3 valued points, or an n-by-3 matrix")

        self._n_points = self._points.shape[0]

        if self._n_points < 2:
            raise ValueError("Path must contain at least 2 points")

        if labels is None:
            self._labels = [f"{self._points[i, 0]:.4g}, {self._points[i, 1]:.4g}, {self._points[i, 2]:.4g}"
                            for i in range(self._n_points)]
        else:
            if len(self._labels) != self._n_points:
                raise ValueError("Expects same number of labels as points")
            self._labels = labels

        self._avoid_endpoints = avoid_endpoints
        self._resolution = resolution

        if scale_by_distance:
            self._section_scalings = []
            for i in range(1, self._n_points):
                self._section_scalings.append(np.sqrt(np.sum((self._points[i] - self._points[i - 1]) ** 2)))
        else:
            self._section_scalings = [1.0 for _ in range(1, self._n_points)]

    def q_points(self):
        """ Get list of q points"""
        output = []
        if self._avoid_endpoints:
            # If we avoid endpoints, we choose points that are close to each of the positions, but not equal

            # Interpolation factor
            f = np.linspace(0, 1, self._resolution).reshape(1, -1)
            f[0,0] = 1e-8
            f[0,-1] = 1-1e-8

            for i in range(1, self._n_points):
                # Linear interpolation
                output.append(self._points[i, :].reshape(-1, 1) * f +
                              self._points[i - 1, :].reshape(-1, 1) * (1 - f))

        else:

            f = np.linspace(0, 1, self._resolution).reshape(1, -1) # Interpolation factor

            # If we don't avoid the endpoints, we need to skip the first point except for the first time
            output.append(self._points[1, :].reshape(-1, 1) * f +
                          self._points[0, :].reshape(-1, 1) * (1 - f))

            for i in range(2, self._n_points):
                # Linear interpolation avoiding skipping the first points
                output.append(self._points[i, :].reshape(-1, 1) * f[0, 1:] +
                              self._points[i - 1, :].reshape(-1, 1) * (1 - f[0, 1:]))

        # Output should be
        #   resolution*(n_points-1), if avoid_endpoints
        #   (resolution-1)*(n_points-1) + 1 otherwise

        return np.concatenate(output, axis=1).T

    def x_values(self):
        """ x values for plotting """
        output = []
        f = np.linspace(0, 1, self._resolution)

        output.append(f * self._section_scalings[0])

        if self._avoid_endpoints:
            for scaling in self._section_scalings[1:]:
                output.append(f*scaling + output[-1][-1])

        else:
            for scaling in self._section_scalings[1:]:
                output.append(f[1:]*scaling + output[-1][-1])

        return np.concatenate(output).T

    def x_ticks(self):
        """ x positions of ticks used to mark positions on a plot"""
        # we want the cumulative sum including a leading zero
        return np.cumsum([0] + self._section_scalings)

    def x_tick_labels(self):
        """ x-axis tick labels corresponding to tick positions """
        return self._labels


if __name__ == "__main__":

    print("Path excluding endpoints, not scaling axes by distance")
    path = Path([[0,0,0],[0,0,1],[1,1,1]], resolution=5)

    print(path.q_points())
    print(path.x_values())
    print(path.x_ticks())
    print(path.x_tick_labels())

    print("Path including endpoints, scaling axes by distance")
    path = Path([[0,0,0],[0,0,1],[1,1,1]],
                resolution=5,
                avoid_endpoints=False,
                scale_by_distance=True)

    print(path.q_points())
    print(path.x_values())
    print(path.x_ticks())
    print(path.x_tick_labels())

