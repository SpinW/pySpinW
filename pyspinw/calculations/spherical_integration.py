""" Support for spherical integration """

from enum import Enum

import numpy as np

from pyspinw.calculations.fibonacci import Fibonacci
from pyspinw.calculations.geodesic import Geodesic

class SphericalPointGenerator:
    """ Base classes for methods of generating points on spheres"""

    method_name = "<base class>"

    def __init__(self, n_points_minimum: int, *args, **kwargs):
        self._actual_n_points = 0
        self._requested_n_points = n_points_minimum
        self._points = np.ndarray((0, 3))
        self._weights = np.ndarray((0, ))

    @property
    def actual_n_points(self) -> int:
        """ Actual number of points produced, always at least as many as requested"""
        return self._actual_n_points

    @property
    def requested_n_points(self) -> int:
        """ Number of points requested """
        return self._requested_n_points

    @property
    def points(self) -> np.ndarray:
        """ 3D Cartesian coordinates of points on sphere"""
        return self._points

    @property
    def weights(self) -> np.ndarray:
        """ weights (dOmega) associated with each point, for integration"""
        return self._weights

    def show_points_lambert(self, new_fig=True, do_show=True):
        """ Show the sample points in the Lambert equal area projection """
        import matplotlib.pyplot as plt

        if new_fig:
            plt.figure("Sample points (Lambert Equal Area) - " + self.method_name)

        # Draw equator and pole

        angles = np.linspace(0, 2*np.pi, 201)
        x = np.sin(angles)
        y = np.cos(angles)

        plt.plot(np.sqrt(2) * x, np.sqrt(2) * y, 'k')
        plt.plot(2 * x, 2 * y, 'k')

        # Plot points

        points = self.points

        z_zero = points[:, 2] == -1

        x = points[:,0] * np.sqrt(2 / (points[:, 2]+1))
        y = points[:,1] * np.sqrt(2 / (points[:, 2]+1))

        x[z_zero] = 0.0
        y[z_zero] = 0.0

        plt.scatter(x, y, color='k')

        if do_show:
            plt.show()



    def show_points_3d(self, new_fig=True, do_show=True):
        """ Show sample points in 3D"""
        import matplotlib.pyplot as plt

        if new_fig:
            fig = plt.figure("Sample points - " + self.method_name)
        else:
            fig = plt.gcf()

        ax = fig.add_subplot(projection='3d')

        points = self.points
        ax.scatter(points[:, 0], points[:, 1], points[:, 2])

        ax.set_xlim3d([-1.1, 1.1])
        ax.set_ylim3d([-1.1, 1.1])
        ax.set_zlim3d([-1.1, 1.1])

        ax.set_box_aspect([1, 1, 1])
        ax.set_proj_type('ortho')

        if do_show:
            plt.show()


class SphericalPointGeneratorType(Enum):
    """ Different kinds of point generators for spherical integration """

    FIBONACCI = "fibonacci"
    RANDOM = "random"
    GEODESIC = "geodesic"

class FibonacciSphericalPointGenerator(SphericalPointGenerator):
    """ Generate points on sphere using the Fibonacci method

    See: J H Hannay and J F Nye 2004 J. Phys. A: Math. Gen. 37 11591
    DOI: 10.1088/0305-4470/37/48/005
    """

    method_name = SphericalPointGeneratorType.FIBONACCI.value.capitalize()

    def __init__(self, n_points_minimum: int, *args, **kwargs):
        super().__init__(n_points_minimum)

        self._points, self._weights = Fibonacci.points_and_weights(n_points_minimum)

        self._actual_n_points = self._weights.shape[0]


class RandomSphericalPointGenerator(SphericalPointGenerator):
    """ Random points on sphere """

    method_name = SphericalPointGeneratorType.RANDOM.value.capitalize()

    def __init__(self, n_points_minimum: int, *args, seed: int | None = None, **kwargs):
        super().__init__(n_points_minimum)

        rng = np.random.default_rng(seed)

        self._actual_n_points = n_points_minimum
        self._points = rng.normal(0,1, size=(n_points_minimum, 3))
        self._points /= np.sqrt(np.sum(self._points**2, axis=1)).reshape((n_points_minimum, 1))
        self._weights = np.ones((n_points_minimum,), dtype=float) * (4 * np.pi / n_points_minimum)


class GeodesicSphericalPointGenerator(SphericalPointGenerator):
    """ Generate points on sphere based on a geodesic geometry"""

    method_name = SphericalPointGeneratorType.GEODESIC.value.capitalize()

    def __init__(self, n_points_minimum: int, *args, **kwargs):
        super().__init__(n_points_minimum)

        divisions = Geodesic.minimal_divisions_for_points(n_points_minimum)
        self._actual_n_points = Geodesic.points_for_division_amount(divisions)

        self._points, self._weights = Geodesic.by_divisions(divisions)

_spherical_point_generator_lookup = {
    cls.method_name.lower(): cls for cls in [
        FibonacciSphericalPointGenerator,
        RandomSphericalPointGenerator,
        GeodesicSphericalPointGenerator]}


def _names_string():
    """ Formatted names """
    all_names = [f"'{name}'" for name in _spherical_point_generator_lookup.keys()]
    return ", ".join(all_names[:-1]) + f" or {all_names[-1]}"


def point_generator(name: SphericalPointGeneratorType | str):
    """ Get a point generator by name / enum value"""
    if isinstance(name, SphericalPointGeneratorType):
        return _spherical_point_generator_lookup[name.value.lower()]

    elif isinstance(name, str):
        try:
            return _spherical_point_generator_lookup[name.lower()]
        except KeyError as ke:
            raise ValueError("Expected `name` to be one of " + _names_string()) from ke

    else:
        raise TypeError("point_generator expected SphericalPointGeneratorType enum or string")

if __name__ == "__main__":

    import matplotlib.pyplot as plt

    for generator in _spherical_point_generator_lookup.values():
        # generator(100).show_points_3d(do_show=False)
        generator(100).show_points_lambert(do_show=False)

    plt.show()
