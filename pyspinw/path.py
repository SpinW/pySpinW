""" Paths through q-space"""
from abc import ABC, abstractmethod

import numpy as np
from numpy._typing import ArrayLike


class Path:
    """ Path through q-space"""

    def __init__(self,
                 points: ArrayLike,
                 n_points_per_segment: int=101,
                 labels: list[str] | None = None,
                 avoid_endpoints=True,
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
        self._n_points_per_segment = n_points_per_segment

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
            f = np.linspace(0, 1, self._n_points_per_segment).reshape(1, -1)
            f[0,0] = 1e-8
            f[0,-1] = 1-1e-8

            for i in range(1, self._n_points):
                # Linear interpolation
                output.append(self._points[i, :].reshape(-1, 1) * f +
                              self._points[i - 1, :].reshape(-1, 1) * (1 - f))

        else:

            f = np.linspace(0, 1, self._n_points_per_segment).reshape(1, -1) # Interpolation factor

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
        f = np.linspace(0, 1, self._n_points_per_segment)

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

    def format_plot(self, plt_or_fig=None):
        """ Apply formatting to a matplotlib plot/figure/axis

        If None, it will import matplotlib.pyplot and work on that
        """
        if plt_or_fig is None:
            import matplotlib.pyplot as plt_or_fig

        if hasattr(plt_or_fig, 'xticks'):
            plt_or_fig.xticks(self.x_ticks(), self.x_tick_labels())
        else:
            plt_or_fig.set_xticks(self.x_ticks(), self.x_tick_labels())

    def __repr__(self):
        path_string = ", ".join([repr(pt) for pt in self._points])
        return f"Path({path_string})"


class ReciprocalSlice:
    """Initialize a 2D slice in reciprocal space.

    Parameters
    ----------
    origin :
        Corner point of the rectangle (3D q-vector)
    a_vec :
        First edge vector (defines the horizontal axis for plotting)
    b_vec :
        Second edge vector (defines the vertical axis for plotting)
    n_a :
        Number of points along a_vec direction
    n_b :
        Number of points along b_vec direction
    labels :
        Optional [a_label, b_label] for plot axes
    """

    def __init__(
        self,
        origin: ArrayLike,
        a_vec: ArrayLike,
        b_vec: ArrayLike,
        n_a: int = 101,
        n_b: int = 101,
        labels: list[str] | None = None,
    ):

        self.origin = np.array(origin)
        self.a_vec = np.array(a_vec)
        self.b_vec = np.array(b_vec)
        self.n_a = np.array(n_a)
        self.n_b = np.array(n_b)

        # Validate shapes
        if self.origin.shape != (3,) or self.a_vec.shape != (3,) or self.b_vec.shape != (3,):
            raise ValueError("origin, a_vec and b_vec must be 3-element vector")

        if labels is None:
            self.labels = ["a (rlu)", "b (rlu)"]
        else:
            if len(labels) != 2:
                raise ValueError("labels must be a list of two strings")
            self.labels = labels

    def q_points(self) -> np.ndarray:
        """Get all q-points in the slice as an (N, 3) array.

        Returns an array of shape (n_a * n_b, 3) with q-points ordered
        """
        # Create parameters s and t from 0 to 1
        s = np.linspace(0, 1, self.n_a)  # shape = (n_a,)
        t = np.linspace(0, 1, self.n_b)  # shape = (n_b,)

        # Create 2D grid of parameters
        s_grid, t_grid = np.meshgrid(s, t, indexing="xy")  # both shape = (n_b, n_a)

        # Flatten to 1D for vectorized calculation
        s_flat = s_grid.flatten()  # shape = (n_a * n_b,)
        t_flat = t_grid.flatten()  # shape = (n_a * n_b,)

        # Calculate q = origin + s*a_vec + t*b_vec (broadcasting)
        q_points = self.origin + s_flat[:, np.newaxis] * self.a_vec + t_flat[:, np.newaxis] * self.b_vec

        return q_points
    
    def grid_shape(self):
          """Shape of the 2D grid for reshaping results: (n_b, n_a).

          This follows matplotlib imshow convention: (rows, columns),
          where n_b is the vertical axis and n_a is the horizontal axis.
          """
          return (self.n_b, self.n_a)
    
    def extent(self) -> list[float]:
          """Extent showing actual projected rlu values along slice axes.

          Returns [a_min, a_max, b_min, b_max] for matplotlib imshow.
          Works correctly for axis-aligned slices like (h,k,0).
          """
          # Find which Cartesian component is dominant in each direction
          a_comp = np.argmax(np.abs(self.a_vec))
          b_comp = np.argmax(np.abs(self.b_vec))

          # Calculate actual min/max values along those projections
          a_min = self.origin[a_comp]
          a_max = self.origin[a_comp] + self.a_vec[a_comp]
          b_min = self.origin[b_comp]
          b_max = self.origin[b_comp] + self.b_vec[b_comp]

          return [a_min, a_max, b_min, b_max]

    def format_plot(self, plt_or_fig=None):
        """Apply x and y labels to a matplotlib plot/figure/axis

        If None, it will import matplotlib.pyplot and work on that
        """
        if plt_or_fig is None:
            import matplotlib.pyplot as plt_or_fig

        if hasattr(plt_or_fig, "xlabel"):
            plt_or_fig.xlabel(self.labels[0])
            plt_or_fig.ylabel(self.labels[1])
        else:
            plt_or_fig.set_xlabel(self.labels[0])
            plt_or_fig.set_ylabel(self.labels[1])

    def __repr__(self):
        return (
            f"ReciprocalSlice(origin={self.origin}, "
            f"a_vec={self.a_vec}, b_vec={self.b_vec}, "
            f"n_a={self.n_a}, n_b={self.n_b})"
        )


class Path1DBase(ABC):
    """ Base class for 1D paths """

    def __init__(self):
        self.n_points = None

    @abstractmethod
    def q_values(self):
        """ Get the q values for this path"""

class Path1D(Path1DBase):
    """ 1D Path, i.e. just values in absolute q """

    def __init__(self,
                 q_min: float = 0.0,
                 q_max: float = 1.0,
                 avoid_endpoints=True,
                 n_points: int = 101):

        super().__init__()

        self.q_min = q_min
        self.q_max = q_max
        self.n_points = n_points
        self.avoid_endpoints = avoid_endpoints

    def q_values(self):
        """ Get q magnitudes """
        base = np.linspace(0, 1, self.n_points)

        if self.avoid_endpoints:
            base[0] = 1e-8
            base[-1] = 1 - 1e-8

        return self.q_min + (self.q_max - self.q_min) * base


class EmpiricalPath1D(Path1DBase):
    """ Path based on data specifying each q value, rather min, max, n"""

    def __init__(self, q_values: ArrayLike):

        super().__init__()

        self._q_values = np.array(q_values)

        if len(self._q_values.shape) != 1:
            raise ValueError("Expected q_values to be a 1D array")

    def q_values(self):
        """ Get the q values"""
        return self._q_values


if __name__ == "__main__":

    print("Path excluding endpoints, not scaling axes by distance")
    path = Path([[0,0,0],[0,0,1],[1,1,1]], n_points_per_segment=5)

    print(path.q_points())
    print(path.x_values())
    print(path.x_ticks())
    print(path.x_tick_labels())

    print("Path including endpoints, scaling axes by distance")
    path = Path([[0,0,0],[0,0,1],[1,1,1]],
                n_points_per_segment=5,
                avoid_endpoints=False,
                scale_by_distance=True)

    print(path.q_points())
    print(path.x_values())
    print(path.x_ticks())
    print(path.x_tick_labels())

