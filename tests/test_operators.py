import numpy as np
import pytest

import spglib

from pyspinw.symmetry.operations import MagneticOperation, SpaceOperation, Operation

r_z = np.array([
    [0, -1,  0],
    [1,  0,  0],
    [0,  0,  1]
])

r_y = np.array([
    [0,  0, -1],
    [0,  1,  0],
    [1,  0,  0]
])

r_x = np.array([
    [1,  0,  0],
    [0,  0, -1],
    [0,  1,  0]
])

f_x = np.array([
    [-1, 0, 0],
    [ 0, 1, 0],
    [ 0, 0, 1]
])


f_y = np.array([
    [1,  0,  0],
    [0, -1,  0],
    [0,  0,  1]
])


f_z = np.array([
    [1,  0,  0],
    [0,  1,  0],
    [0,  0, -1]
])



def random_rotoreflection(rng: np.random.Generator):
    r = np.eye(3)
    for f in [f_x, f_y, f_z]:
        if rng.random() < 0.5:
            r @= f

    for rotation in [r_x, r_y, r_z]:
        for i in range(rng.integers(0, 4)):
            r @= rotation

    for f in [f_x, f_y, f_z]:
        if rng.random() < 0.5:
            r @= f

    return r


rng = np.random.default_rng(8746)
test_points = [rng.random((10, 3)) for i in range(5)]
test_spins = [rng.random((10, 3)) for i in range(5)]
test_points_and_spins = [tup for tup in zip(test_points, test_spins)]

random_magnetic_operation_info = [(random_rotoreflection(rng),
                                   rng.random((3,)),
                                   rng.integers(0, 2) * 2 - 1)
                                  for i in range(10)]

random_space_operation_info = [(random_rotoreflection(rng),
                                   rng.random((3,)))
                                  for i in range(10)]

@pytest.mark.parametrize("identity", [
    SpaceOperation.from_numpy(np.eye(3), np.zeros(3), name="Id"),
    MagneticOperation.from_numpy(np.eye(3), np.zeros(3), 1, "Id")
])
@pytest.mark.parametrize("points_and_spins", test_points_and_spins)
def test_identity(identity: Operation, points_and_spins: tuple[np.ndarray, np.ndarray]):
    """ Check that the identity operation does nothing """

    points, spins = points_and_spins

    comparison_points = identity.transform_positions(points)
    assert np.allclose(points, comparison_points)

    comparison_points, comparison_spins = identity.transform_positions_and_spins(points, spins)
    assert np.allclose(points, comparison_points)
    assert np.allclose(spins, comparison_spins)


@pytest.mark.parametrize("points_and_spins", test_points_and_spins)
def test_time_reversal_only(points_and_spins):
    """ Check that time-reversal only inverts moments"""

    points, spins = points_and_spins

    generator = MagneticOperation.from_numpy(np.eye(3), np.zeros(3), -1, "Id'")
    comparison_points, comparison_spins = generator.transform_positions_and_spins(points, spins)

    assert np.allclose(points, comparison_points)
    assert np.allclose(spins, -comparison_spins)

@pytest.mark.parametrize("cls, params, time_reverse", [
    (SpaceOperation, {}, 1),
    (MagneticOperation, {"time_reversal": -1}, -1),
    (MagneticOperation, {"time_reversal": 1}, 1)])
@pytest.mark.parametrize("k", [0.2, 0.4, 0.8])
@pytest.mark.parametrize("points_and_spins", test_points_and_spins)
def test_translation(cls, params, points_and_spins, k, time_reverse):
    """ Translation only"""

    points, spins = points_and_spins

    operation = cls.from_numpy(np.eye(3), np.zeros(3)+k, **params)

    point_comparison = operation.transform_positions(points)
    assert np.allclose((points + k) % 1, point_comparison)

    point_comparison, spin_comparison = operation.transform_positions_and_spins(points, spins)
    assert np.allclose((points + k) % 1, point_comparison)
    assert np.allclose(time_reverse * spins, spin_comparison)



@pytest.mark.parametrize("points, spins", test_points_and_spins)
@pytest.mark.parametrize("rotation, translation, reversal", [
    (np.eye(3)[(0,2,1), :], np.zeros(3), 1), # Swap y and z
    (np.eye(3), np.zeros(3)+0.5, 1), # Translate half a cell
    (np.eye(3), np.zeros(3), -1)  ]) # Time reverse
def test_z2_magnetic(points, spins, rotation, translation, reversal):
    """ Check period 2 generators have f(f(x)) = x"""
    g = MagneticOperation.from_numpy(rotation, translation, reversal)

    comparison_points = g.transform_positions(g.transform_positions(points))
    assert np.allclose(points, comparison_points)

    comparison_points, comparison_spins = g.transform_positions_and_spins(
                                            *g.transform_positions_and_spins(
                                                points, spins))
    assert np.allclose(points, comparison_points)
    assert np.allclose(spins, comparison_spins)


@pytest.mark.parametrize("points, spins", test_points_and_spins)
@pytest.mark.parametrize("rotation, translation", [
    (np.eye(3)[(0,2,1), :], np.zeros(3)), # Swap y and z
    (np.eye(3), np.zeros(3)+0.5)]) # Translate half a cell
def test_z2_space(points, spins, rotation, translation):
    """ Check period 2 generators have f(f(x)) = x"""
    g = SpaceOperation.from_numpy(rotation, translation)
    comparison_points = g.transform_positions(g.transform_positions(points))

    assert np.allclose(points, comparison_points)

    comparison_points, comparison_spins = g.transform_positions_and_spins(
        *g.transform_positions_and_spins(
            points, spins))

    assert np.allclose(points, comparison_points)
    assert np.allclose(spins, comparison_spins)




@pytest.mark.parametrize("points, spins", test_points_and_spins)
@pytest.mark.parametrize("translation, rotation, reversal", [
    (np.eye(3)[(0,2,1), :], np.zeros(3), 1), # Swap y and z
    (np.eye(3), np.zeros(3)+0.5, 1), # Translate half a cell
    (np.eye(3), np.zeros(3), -1)  ]) # Time reverse
def test_z2_composition_magnetic(points, spins, translation, rotation, reversal):
    """ Check some generators that have f(f(x)) = x, using composition"""
    g = MagneticOperation.from_numpy(translation, rotation, reversal)

    comparison_points = g.and_then(g).transform_positions(points)
    assert np.allclose(points, comparison_points)

    comparison_points, comparison_spins = g.and_then(g).transform_positions_and_spins(points, spins)
    assert np.allclose(points, comparison_points)
    assert np.allclose(spins, comparison_spins)


@pytest.mark.parametrize("points, spins", test_points_and_spins)
@pytest.mark.parametrize("translation, rotation", [
    (np.eye(3)[(0,2,1), :], np.zeros(3)), # Swap y and z
    (np.eye(3), np.zeros(3)+0.5)]) # Translate half a cell
def test_z2_composition_space(points, spins, translation, rotation):
    """ Check some generators that have f(f(x)) = x, using composition"""
    g = SpaceOperation.from_numpy(translation, rotation)

    comparison_points = g.and_then(g).transform_positions(points)
    assert np.allclose(points, comparison_points)

    comparison_points, comparison_spins = g.and_then(g).transform_positions_and_spins(points, spins)
    assert np.allclose(points, comparison_points)
    assert np.allclose(spins, comparison_spins)


@pytest.mark.parametrize("points, spins", test_points_and_spins)
@pytest.mark.parametrize("g1", random_magnetic_operation_info)
@pytest.mark.parametrize("g2", random_magnetic_operation_info)
def test_magnetic_composition_no_translation_or_reversal(points, spins, g1, g2):
    generator_1 = MagneticOperation.from_numpy(g1[0], np.zeros(3), 1, name="Generator 1")
    generator_2 = MagneticOperation.from_numpy(g2[0], np.zeros(3), 1, name="Generator 2")

    transformed_points_sequential, transformed_spins_sequential = \
        generator_2.transform_positions_and_spins(
            *generator_1.transform_positions_and_spins(points, spins))

    transformed_points_compose, transformed_spins_compose = \
        generator_1.and_then(generator_2).transform_positions_and_spins(points, spins)

    assert np.allclose(transformed_points_compose, transformed_points_sequential)
    assert np.allclose(transformed_spins_compose, transformed_spins_sequential)


@pytest.mark.parametrize("points, spins", test_points_and_spins)
@pytest.mark.parametrize("g1", random_magnetic_operation_info)
@pytest.mark.parametrize("g2", random_magnetic_operation_info)
def test_space_composition_no_translation_or_reversal(points, spins, g1, g2):
    generator_1 = SpaceOperation.from_numpy(g1[0], np.zeros(3), name="Generator 1")
    generator_2 = SpaceOperation.from_numpy(g2[0], np.zeros(3), name="Generator 2")

    transformed_points_sequential, transformed_spins_sequential = \
        generator_2.transform_positions_and_spins(
            *generator_1.transform_positions_and_spins(points, spins))

    transformed_points_compose, transformed_spins_compose = \
        generator_1.and_then(generator_2).transform_positions_and_spins(points, spins)

    assert np.allclose(transformed_points_compose, transformed_points_sequential)
    assert np.allclose(transformed_spins_compose, transformed_spins_sequential)


@pytest.mark.parametrize("points, spins", test_points_and_spins)
@pytest.mark.parametrize("g1", random_magnetic_operation_info)
@pytest.mark.parametrize("g2", random_magnetic_operation_info)
def test_magnetic_composition_full(points, spins, g1, g2):
    generator_1 = MagneticOperation.from_numpy(g1[0], g1[1], g1[2], name="Generator 1")
    generator_2 = MagneticOperation.from_numpy(g2[0], g2[1], g2[2], name="Generator 2")

    transformed_points_sequential, transformed_spins_sequential = \
        generator_2.transform_positions_and_spins(
            *generator_1.transform_positions_and_spins(points, spins))

    transformed_points_compose, transformed_spins_compose = \
        generator_1.and_then(generator_2).transform_positions_and_spins(points, spins)

    assert np.allclose(transformed_points_compose, transformed_points_sequential)
    assert np.allclose(transformed_spins_compose, transformed_spins_sequential)


@pytest.mark.parametrize("points, spins", test_points_and_spins)
@pytest.mark.parametrize("g1", random_magnetic_operation_info)
@pytest.mark.parametrize("g2", random_magnetic_operation_info)
def test_space_composition_full(points, spins, g1, g2):
    generator_1 = SpaceOperation.from_numpy(g1[0], g1[1], name="Generator 1")
    generator_2 = SpaceOperation.from_numpy(g2[0], g2[1], name="Generator 2")

    transformed_points_sequential, transformed_spins_sequential = \
        generator_2.transform_positions_and_spins(
            *generator_1.transform_positions_and_spins(points, spins))

    transformed_points_compose, transformed_spins_compose = \
        generator_1.and_then(generator_2).transform_positions_and_spins(points, spins)

    assert np.allclose(transformed_points_compose, transformed_points_sequential)
    assert np.allclose(transformed_spins_compose, transformed_spins_sequential)

@pytest.mark.parametrize("g", random_magnetic_operation_info)
def test_magnetic_operation_equality_but_not_identity(g):

    generator_1 = MagneticOperation.from_numpy(*g, name="Generator 1")
    generator_2 = MagneticOperation.from_numpy(*g, name="Generator 2")

    assert generator_1 is not generator_2
    assert generator_1 == generator_2


@pytest.mark.parametrize("g", random_space_operation_info)
def test_space_operation_equality_but_not_identity(g):

    generator_1 = SpaceOperation.from_numpy(*g, name="Generator 1")
    generator_2 = SpaceOperation.from_numpy(*g, name="Generator 2")

    assert generator_1 is not generator_2
    assert generator_1 == generator_2


@pytest.mark.parametrize("index_1", range(len(random_magnetic_operation_info)))
@pytest.mark.parametrize("delta", range(len(random_magnetic_operation_info) - 1))
def test_magnetic_operation_not_equal(index_1: int, delta: int):
    """ Check non-equal pairs are not identified as equal"""
    index_2 = (index_1 + 1 + delta) % len(random_magnetic_operation_info)

    g1 = random_magnetic_operation_info[index_1]
    g2 = random_magnetic_operation_info[index_2]

    generator_1 = MagneticOperation.from_numpy(*g1, name="Generator 1")
    generator_2 = MagneticOperation.from_numpy(*g2, name="Generator 2")

    assert generator_1 is not generator_2
    assert generator_1 != generator_2


@pytest.mark.parametrize("index_1", range(len(random_space_operation_info)))
@pytest.mark.parametrize("delta", range(len(random_space_operation_info) - 1))
def test_space_operation_not_equal(index_1: int, delta: int):
    """ Check non-equal pairs are not identified as equal"""
    index_2 = (index_1 + 1 + delta) % len(random_space_operation_info)

    g1 = random_space_operation_info[index_1]
    g2 = random_space_operation_info[index_2]

    generator_1 = SpaceOperation.from_numpy(*g1, name="Generator 1")
    generator_2 = SpaceOperation.from_numpy(*g2, name="Generator 2")

    assert generator_1 is not generator_2
    assert generator_1 != generator_2
