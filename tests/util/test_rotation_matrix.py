import pytest
import numpy as np

from pyspinw.util import rotation_matrix

@pytest.mark.parametrize("start,axis,expected", [
    [[1, 1, 1], [ 1, 0, 0], [1, 1,-1]],
    [[1, 1, 1], [-1, 0, 0], [1,-1, 1]],
    [[1, 0, 0], [ 1, 0, 0], [1, 0, 0]],
    [[0, 1, 0], [ 1, 0, 0], [0, 0,-1]],
    [[0, 1, 0], [-1, 0, 0], [0, 0, 1]],
    [[0, 1, 0], [ 0, 0, 1], [ 1, 0, 0]],
    [[0, 1, 0], [ 0, 0,-1], [-1, 0, 0]],
    ])
def test_90_rotation_in_cartesian(start, axis, expected):
    """ Check that 90 deg rotations work """
    assert start @ rotation_matrix(np.pi/2, axis) == pytest.approx(expected)


@pytest.mark.parametrize("angle", [1,2,3])
@pytest.mark.parametrize("scale", [1,2,3])
def test_axis_scaling_invariance(angle, scale):
    """ Check that the length of the axis vector doesn't matter"""
    assert rotation_matrix(angle, [scale, scale, scale]) == pytest.approx(rotation_matrix(angle, [1, 1, 1]))

def test_bad_axis_raises_error():
    """ Check that axis=[0,0,0] raises an exception """
    with pytest.raises(ValueError):
        rotation_matrix(1, [0,0,0])


vectors = [
    [1,0,0],
    [1,1,1],
    [1,2,3],
    [-1,1,-1],
    [-10,7,1]]

@pytest.mark.parametrize("angle", [0.1, 0.5, 1.0])
@pytest.mark.parametrize("axis_vector", vectors)
@pytest.mark.parametrize("vector", vectors)
def test_rotation_preserves_length(angle, axis_vector, vector):
    """ Check that the rotation matrices preserve length"""
    axis_vector = np.array(axis_vector)
    vector = np.array(vector)

    if np.all(axis_vector == vector):
        axis = axis_vector
    else:
        axis = np.cross(vector, axis_vector)

    rotated = vector @ rotation_matrix(angle, axis)

    assert np.sum(np.array(vector**2)) == pytest.approx(np.sum(rotated**2))


@pytest.mark.parametrize("angle", [0.1, 0.5, 1.0])
@pytest.mark.parametrize("axis_vector", vectors)
@pytest.mark.parametrize("vector", vectors)
def test_rotation_amount(angle, axis_vector, vector):
    """ Check that the angle rotated by the rotation matrix is the correct angle"""
    axis_vector = np.array(axis_vector)
    vector = np.array(vector)

    # Get an axis perpendicular to the input vector
    if np.all(axis_vector == vector):
        axis = np.cross(vector, [0, 0, 1])
    else:
        axis = np.cross(vector, axis_vector)

    rotated = vector @ rotation_matrix(angle, axis)

    mag = np.sqrt(np.sum(rotated**2)*np.sum(vector**2))

    assert np.dot(rotated, vector) == pytest.approx(mag*np.abs(np.cos(angle)))

