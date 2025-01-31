import numpy as np
import pytest

from pyspinw.util.group_generators import Generator

def random_rotoreflection(rng: np.random.Generator):
    # Can only have -1,0,1 in each location
    return rng.integers(-1, 2, size=(3, 3))

rng = np.random.default_rng(8746)
test_points = [rng.random((2, 6)) for i in range(5)]
random_generator_info = [(random_rotoreflection(rng),
                          rng.random((3,)),
                          rng.integers(0, 2) * 2 - 1)
                         for i in range(5)]

@pytest.mark.parametrize("points", test_points)
def test_identity(points):
    """ Check that the identity generator does nothing """
    generator = Generator(np.eye(3), np.zeros(3), 1, "Id")
    comparison = generator(points)

    assert np.all(np.abs(points - comparison) < 1e-10)


@pytest.mark.parametrize("points", test_points)
def test_time_reversal_only(points):
    """ Check that time-reversal only inverts moments"""
    generator = Generator(np.eye(3), np.zeros(3), -1, "Id")
    comparison = generator(points)

    assert np.all(np.abs(points[:, :3] - comparison[:, :3]) < 1e-10)
    assert np.all(np.abs(points[:, 3:] + comparison[:, 3:]) < 1e-10)


@pytest.mark.parametrize("k", [0.2, 0.4, 0.8])
@pytest.mark.parametrize("points", test_points)
def test_translation(points, k):
    """ Translation only"""
    generator = Generator(np.eye(3), np.zeros(3) + k, -1, "Id")
    comparison = generator(points)

    assert np.all(np.abs((points[:, :3] + k)%1 - comparison[:, :3]) < 1e-10)


@pytest.mark.parametrize("points", test_points)
@pytest.mark.parametrize("generator", [
    (np.eye(3)[(0,2,1), :], np.zeros(3), 1), # Swap y and z
    (np.eye(3), np.zeros(3)+0.5, 1), # Translate half a cell
    (np.eye(3), np.zeros(3), -1)  ]) # Time reverse
def test_z2(points, generator):
    """ Check some generators that have f(f(x)) = x"""
    g = Generator(generator[0], generator[1], generator[2], "Id")
    comparison = g(g(points))

    assert np.all(np.abs(points - comparison) < 1e-10)



@pytest.mark.parametrize("points", test_points)
@pytest.mark.parametrize("generator", [
    (np.eye(3)[(0,2,1), :], np.zeros(3), 1), # Swap y and z
    (np.eye(3), np.zeros(3)+0.5, 1), # Translate half a cell
    (np.eye(3), np.zeros(3), -1)  ]) # Time reverse
def test_z2_composition(points, generator):
    """ Check some generators that have f(f(x)) = x, using composition"""
    g = Generator(generator[0], generator[1], generator[2], "Id")
    comparison = g.and_then(g)(points)

    assert np.all(np.abs(points - comparison) < 1e-10)


@pytest.mark.parametrize("points", test_points)
@pytest.mark.parametrize("g1", random_generator_info)
@pytest.mark.parametrize("g2", random_generator_info)
def test_composition_no_translation_or_reversal(points, g1, g2):
    generator_1 = Generator(g1[0], np.zeros(3), 1, name="Generator 1")
    generator_2 = Generator(g2[0], np.zeros(3), 1, name="Generator 2")

    transformed_points_sequential = generator_1(generator_2(points))

    # transformed_points_compose = generator_1.and_then(generator_2)(points)
    transformed_points_compose = generator_2.and_then(generator_1)(points)

    assert np.all(np.abs(transformed_points_compose - transformed_points_sequential) < 1e-10)



@pytest.mark.parametrize("points", test_points)
@pytest.mark.parametrize("g1", random_generator_info)
@pytest.mark.parametrize("g2", random_generator_info)
def test_composition_full(points, g1, g2):
    generator_1 = Generator(*g1, name="Generator 1")
    generator_2 = Generator(*g2, name="Generator 2")

    transformed_points_sequential = generator_2(generator_1(points))
    transformed_points_compose = generator_1.and_then(generator_2)(points)

    assert np.all(np.abs(transformed_points_compose - transformed_points_sequential) < 1e-10)

