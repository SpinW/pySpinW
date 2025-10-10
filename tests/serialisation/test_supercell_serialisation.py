""" Tests for the serialisation of supercells """
import numpy as np

from pyspinw.symmetry.supercell import TrivialSupercell, Supercell, \
    CommensuratePropagationVector, SummationSupercell, TransformationSupercell, RotationTransform


def test_trivial_supercell_serialisation():
    """ Check that the trivial supercell works """
    supercell = TrivialSupercell(scaling=(1,2,3))
    json = supercell.serialise()
    deserialised = Supercell.deserialise(json)

    assert isinstance(deserialised, TrivialSupercell)

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


def test_rotation_supercell_serialisation():
    """ Test that rotation supercells serialise correctly"""
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

