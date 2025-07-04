import numpy as np
import pytest

import spglib

from pyspinw.symmetry.operations import MagneticOperation

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
test_points = [rng.random((10, 6)) for i in range(5)]

random_generator_info = [(random_rotoreflection(rng),
                          rng.random((3,)),
                          rng.integers(0, 2) * 2 - 1)
                         for i in range(10)]

@pytest.mark.parametrize("points", test_points)
def test_identity(points):
    """ Check that the identity generator does nothing """
    generator = MagneticOperation.from_numpy(np.eye(3), np.zeros(3), 1, "Id")
    comparison = generator(points)

    assert np.all(np.abs(points - comparison) < 1e-10)


@pytest.mark.parametrize("points", test_points)
def test_time_reversal_only(points):
    """ Check that time-reversal only inverts moments"""
    generator = MagneticOperation.from_numpy(np.eye(3), np.zeros(3), -1, "Id'")
    comparison = generator(points)

    assert np.all(np.abs(points[:, :3] - comparison[:, :3]) < 1e-10)
    assert np.all(np.abs(points[:, 3:] + comparison[:, 3:]) < 1e-10)


@pytest.mark.parametrize("k", [0.2, 0.4, 0.8])
@pytest.mark.parametrize("points", test_points)
def test_translation(points, k):
    """ Translation only"""
    generator = MagneticOperation.from_numpy(np.eye(3), np.zeros(3)+k, -1)
    comparison = generator(points)

    assert np.all(np.abs((points[:, :3] + k)%1 - comparison[:, :3]) < 1e-10)


@pytest.mark.parametrize("points", test_points)
@pytest.mark.parametrize("generator", [
    (np.eye(3)[(0,2,1), :], np.zeros(3), 1), # Swap y and z
    (np.eye(3), np.zeros(3)+0.5, 1), # Translate half a cell
    (np.eye(3), np.zeros(3), -1)  ]) # Time reverse
def test_z2(points, generator):
    """ Check some generators that have f(f(x)) = x"""
    g = MagneticOperation.from_numpy(generator[0], generator[1], generator[2], "Id")
    comparison = g(g(points))

    assert np.all(np.abs(points - comparison) < 1e-10)



@pytest.mark.parametrize("points", test_points)
@pytest.mark.parametrize("generator", [
    (np.eye(3)[(0,2,1), :], np.zeros(3), 1), # Swap y and z
    (np.eye(3), np.zeros(3)+0.5, 1), # Translate half a cell
    (np.eye(3), np.zeros(3), -1)  ]) # Time reverse
def test_z2_composition(points, generator):
    """ Check some generators that have f(f(x)) = x, using composition"""
    g = MagneticOperation.from_numpy(generator[0], generator[1], generator[2])
    comparison = g.and_then(g)(points)

    assert np.all(np.abs(points - comparison) < 1e-10)


@pytest.mark.parametrize("points", test_points)
@pytest.mark.parametrize("g1", random_generator_info)
@pytest.mark.parametrize("g2", random_generator_info)
def test_composition_no_translation_or_reversal(points, g1, g2):
    generator_1 = MagneticOperation.from_numpy(g1[0], np.zeros(3), 1, name="Generator 1")
    generator_2 = MagneticOperation.from_numpy(g2[0], np.zeros(3), 1, name="Generator 2")

    transformed_points_sequential = generator_1(generator_2(points))

    # transformed_points_compose = generator_1.and_then(generator_2)(points)
    transformed_points_compose = generator_2.and_then(generator_1)(points)

    assert np.all(np.abs(transformed_points_compose - transformed_points_sequential) < 1e-10)



@pytest.mark.parametrize("points", test_points)
@pytest.mark.parametrize("g1", random_generator_info)
@pytest.mark.parametrize("g2", random_generator_info)
def test_composition_full(points, g1, g2):
    """ Test that g2(g1(x)) == (g1 and_then g2)(x) """
    generator_1 = MagneticOperation.from_numpy(*g1, name="Generator 1")
    generator_2 = MagneticOperation.from_numpy(*g2, name="Generator 2")

    transformed_points_sequential = generator_2(generator_1(points))
    transformed_points_compose = generator_1.and_then(generator_2)(points)

    assert np.all(np.abs(transformed_points_compose - transformed_points_sequential) < 1e-10)

#
# generators_for_testing = _spglib_generators_to_objects(spglib.get_magnetic_symmetry_from_database(1651))
# seeds = [1234, 76423, 2093478, 7973634]
# @pytest.mark.parametrize("seed", seeds)
# def test_generator_sorting(seed: int):
#     rng = np.random.default_rng(seed)
#     random_order = np.arange(len(generators_for_testing))
#     rng.shuffle(random_order)
#
#     initial_copy = generators_for_testing.copy()
#     other_copy = [initial_copy[i] for i in random_order]
#
#     initial_copy.sort()
#     other_copy.sort()
#
#     for a, b in zip(initial_copy, other_copy):
#         assert a == b

@pytest.mark.parametrize("g", random_generator_info)
def test_generator_equality_but_not_identity(g):

    generator_1 = MagneticOperation.from_numpy(*g, name="Generator 1")
    generator_2 = MagneticOperation.from_numpy(*g, name="Generator 2")

    assert generator_1 is not generator_2
    assert generator_1 == generator_2


@pytest.mark.parametrize("index_1", range(len(random_generator_info)))
@pytest.mark.parametrize("delta", range(len(random_generator_info)-1))
def test_generator_not_equal(index_1: int, delta: int):
    """ Check non-equal pairs are not identified as equal"""
    index_2 = (index_1 + 1 + delta) % len(random_generator_info)

    g1 = random_generator_info[index_1]
    g2 = random_generator_info[index_2]

    generator_1 = MagneticOperation.from_numpy(*g1, name="Generator 1")
    generator_2 = MagneticOperation.from_numpy(*g2, name="Generator 2")

    assert generator_1 is not generator_2
    assert generator_1 != generator_2
