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


