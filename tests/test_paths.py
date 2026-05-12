""" Tests for the Path class """

import numpy as np
import pytest

from pyspinw.path import Path, Slice

point_list = [[0,0,0], [0,0,1], [0,1,0], [1,0,0], [1,1,1]]

@pytest.mark.parametrize("n_points", [2,3,4,5])
@pytest.mark.parametrize("avoid_endpoints", [True, False])
def test_sizes_match(n_points, avoid_endpoints):
    """ Check the sizes of the outputs of the different methods are consistent """
    resolution = 101
    path = Path(point_list[:n_points],
                avoid_endpoints=avoid_endpoints,
                n_points_per_segment=resolution)

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


@pytest.mark.parametrize("avoid_endpoints", [False, True])
def test_path_slices(avoid_endpoints):
    critical_points = [[0,0,0], [0,0,1], [0,1,1], [0,1,0]]
    path = Path(critical_points, avoid_endpoints=avoid_endpoints)

    points = path.q_points()

    assert len(critical_points) == len(path.section_slices()) + 1, "There should be n_points-1 slices"

    for s in path.section_slices():
        test_points = points[s, :]

        deltas = test_points[1:, :] - test_points[:-1, :]

        for i in range(1, deltas.shape[0]):
            assert np.allclose(deltas[i-1, :], deltas[i, :]), "Direction should be continuous along sections"


def test_slice_q_points_shape():
    """q_points returns correct shape for different n_a, n_b"""
    s = Slice([0,0,0], [1,0,0], [0,1,0], n_a=3, n_b=5)
    assert s.q_points().shape == (15, 3)


def test_slice_grid_shape():
    """grid_shape returns (n_a, n_b)"""
    s = Slice([0,0,0], [1,0,0], [0,1,0], n_a=3, n_b=5)
    assert s.grid_shape() == (3, 5)


def test_slice_padding():
    """padding extends the q range correctly"""
    s = Slice([0,0,0], [1,0,0], [0,1,0], n_a=5, padding=0.1)
    q = s.q_points()
    assert np.min(q[:,0]) == -0.1
    assert np.max(q[:,0]) == 1.1
