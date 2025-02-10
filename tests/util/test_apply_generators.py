import pytest

import numpy as np
from pyspinw.util.apply_generators import _remove_duplicates_and_order_points

def test_order_unique_helper_removes_duplicates():
    """ Test that the helper function that finds unique entries removes duplicates"""

    test_data = np.array([
        [1,2,3],
        [3,4,5],
        [2,3,4],
        [2,3,4],
        [1,2,3],
        [1,2,3]
    ])

    expected = np.array([
        [1,2,3],
        [2,3,4],
        [3,4,5]
    ])

    unique = _remove_duplicates_and_order_points(test_data)

    assert unique.shape == (3, 3)

    assert unique == pytest.approx(expected, abs=1e-10)

def test_order_unique_helper_gives_well_defined_order():
    """ Test that ordering of points is uniquely defined """
    test_data = np.array([
        [1,1,1],
        [1,0,1],
        [1,0,1],
        [0,1,1],
        [0,0,0], # 5
        [0,1,0],
        [0,1,1],
        [1,1,0],
        [0,0,1],
        [0,1,0]  # 10
    ])

    another_order = [3,2,7,5,9,1,0,4,8,6]

    other_test_data = test_data[another_order, :]

    first_way = _remove_duplicates_and_order_points(test_data)
    second_way = _remove_duplicates_and_order_points(other_test_data)

    assert first_way == pytest.approx(second_way, abs=1e-10)