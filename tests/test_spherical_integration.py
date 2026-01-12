import pytest

import numpy as np

from pyspinw.calculations.spherical_integration import _spherical_point_generator_lookup


point_generator_classes = [cls for cls in _spherical_point_generator_lookup.values()]


@pytest.mark.parametrize("cls", point_generator_classes)
def test_weight_sum(cls):
    """ Check that the weights provided by integration methods add up to 4 pi"""
    generator = cls(n_points_minimum=100)
    assert np.sum(generator.weights) == pytest.approx(4 * np.pi, abs=1e-9), "Weights should add up to 4pi"

@pytest.mark.parametrize("cls", point_generator_classes)
@pytest.mark.parametrize("n_points_minimum", [10, 20, 50, 100, 200, 500, 1000])
def test_points_shape(cls, n_points_minimum):
    """ Check that there are enough points, and that the output is the right shape"""
    generator = cls(n_points_minimum=n_points_minimum)

    assert generator.actual_n_points >= n_points_minimum, "Should be at least the number of points requested"
    assert generator.points.shape == (generator.actual_n_points, 3), "Output shape should be actual_points-by-3"

@pytest.mark.parametrize("cls", point_generator_classes)
def test_points_are_on_sphere(cls):
    """ Check that there are enough points, and that the output is the right shape"""
    generator = cls(n_points_minimum=1000)

    points = generator.points

    assert np.all(np.isclose(np.sum(points**2, axis=1), np.ones((generator.actual_n_points, )), 1e-9))



def f0(x,y,z):
    """ Integrates to 4 pi"""
    return np.ones_like(x)

def f1(x,y,z):
    """ Integrates to zero """
    return x*y*z

def f2(x,y,z):
    """ Integrates to pi/35 """
    return (x*y*z)**2

def f3(x,y,z):
    """ Integrates to 41 pi / 15"""
    return x**4 + y**4 + z**4

test_functions = [
    (f0, 4 * np.pi),
    (f1, 0),
    (f2, 4 * np.pi / 105),
    (f3, 12 * np.pi / 5)]

@pytest.mark.parametrize("function, integral", test_functions)
@pytest.mark.parametrize("cls", point_generator_classes)
def test_integrate_functions(cls, function, integral):
    """ Check integration with known function, f(x,y,z)= (x y z)^2, integral should be pi/35"""

    if cls.method_name == "Random":
        n = 100_000
        tol = 1e-2
    else:
        n = 1000
        tol = 1e-4

    generator = cls(n_points_minimum=n)

    points = generator.points
    sample_values = function(points[:, 0], points[:, 1], points[:, 2])

    numerical_integral = np.sum(sample_values * generator.weights)

    assert numerical_integral == pytest.approx(integral, abs=tol), \
        "Numerical integral should be approximately the analytical integral"
