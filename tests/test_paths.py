""" Tests for the Path class """

import pytest

from pyspinw.path import Path

point_list = [[0,0,0], [0,0,1], [0,1,0], [1,0,0], [1,1,1]]

@pytest.mark.parametrize("n_points", [2,3,4,5])
@pytest.mark.parametrize("avoid_endpoints", [True, False])
def test_sizes_match(n_points, avoid_endpoints):
    """ Check the sizes of the outputs of the different methods are consistent """
    resolution = 101
    path = Path(point_list[:n_points],
                avoid_endpoints=avoid_endpoints,
                resolution=resolution)

    # Check the graph x and q value sizes match
    assert len(path.q_points()) == len(path.x_values())

    # Check the lengths are what we would predict
    if avoid_endpoints:
        assert len(path.x_values()) == (n_points-1)*resolution

    else:
        assert len(path.x_values()) == (n_points-1)*(resolution-1) + 1

    # check the labelling stuff is right
    assert len(path.x_ticks()) == len(path.x_tick_labels())
    assert len(path.x_ticks()) == n_points


