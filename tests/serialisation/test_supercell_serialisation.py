""" Tests for the serialisation of supercells """
import numpy as np

from pyspinw.symmetry.supercell import TiledSupercell, Supercell, \
    CommensuratePropagationVector, SummationSupercell, TransformationSupercell, RotationTransform, \
    DirectSupercell, RotationSupercell, PropagationVector


def test_trivial_supercell_serialisation():
    """ Check that the trivial supercell works """
    supercell = TiledSupercell(scaling=(1, 2, 3))
    json = supercell.serialise()
    deserialised = Supercell.deserialise(json)

    assert isinstance(deserialised, TiledSupercell)

    assert deserialised._scaling == supercell._scaling

def test_summation_supercell_serialisation():
    """ Check that summation supercells serialise correctly """
    vectors = [CommensuratePropagationVector(0, 0, 1/2),
               CommensuratePropagationVector(1/3, 1/3, 1/3)]
    supercell = SummationSupercell(vectors, scaling=(1,3,5))

    json = supercell.serialise()
    deserialised = Supercell.deserialise(json)

    assert isinstance(deserialised, SummationSupercell)

    assert deserialised._scaling == supercell._scaling

    assert all([preserialised_vector == deserialised_vector
               for preserialised_vector, deserialised_vector
                in zip(supercell._propagation_vectors, deserialised._propagation_vectors)])


def test_transformation_supercell_serialisation():
    """ Test that transformation supercells serialise correctly"""
    input_data = [(CommensuratePropagationVector(0, 0, 1 / 2), RotationTransform([0,1,0])),
                (CommensuratePropagationVector(1 / 3, 1 / 3, 1 / 3), RotationTransform([1,0,0]))]
    supercell = TransformationSupercell(input_data, scaling=(1, 3, 5))

    json = supercell.serialise()
    deserialised = Supercell.deserialise(json)

    assert isinstance(deserialised, TransformationSupercell)

    assert deserialised._scaling == supercell._scaling

    assert all([preserialised_vector == deserialised_vector
                for preserialised_vector, deserialised_vector
                in zip(supercell._propagation_vectors, deserialised._propagation_vectors)])

    for (_, unserialised), (_, deserialised) in zip(supercell._transforms, deserialised._transforms):
        assert isinstance(deserialised, RotationTransform)
        assert isinstance(unserialised, RotationTransform) # Should be true by this test
        assert np.all(np.abs(deserialised._axis - unserialised._axis) < 1e-10)

def test_direct_supercell_serialisation():
    """ Check that the direct supercell serialises correctly """
    supercell = DirectSupercell(5,6,7, scaling=(10,11,12))

    json = supercell.serialise()
    deserialised = Supercell.deserialise(json)

    assert isinstance(deserialised, DirectSupercell)

    assert supercell.a == deserialised.a
    assert supercell.b == deserialised.b
    assert supercell.c == deserialised.c

    assert supercell.scaling == deserialised.scaling


def test_rotation_supercell_serialisation():
    """ Check that rotation supercells serialise correctly """

    supercell = RotationSupercell([2, 3, 4],
                                  PropagationVector(1 / np.sqrt(103), 1 / np.sqrt(102), 1 / np.sqrt(101),
                                                    phase=np.pi/3))

    json = supercell.serialise()

    print(json)

    deserialised = Supercell.deserialise(json)

    assert isinstance(deserialised, RotationSupercell)

    assert supercell.propagation_vector.i == deserialised.propagation_vector.i
    assert supercell.propagation_vector.j == deserialised.propagation_vector.j
    assert supercell.propagation_vector.k == deserialised.propagation_vector.k
    assert supercell.propagation_vector.phase == deserialised.propagation_vector.phase

    assert np.allclose(supercell.perpendicular, deserialised.perpendicular)