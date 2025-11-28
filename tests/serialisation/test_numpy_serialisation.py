import pytest
import numpy as np

from pyspinw.serialisation import numpy_serialise, numpy_deserialise, SPWSerialisationError

int_examples = [
    np.arange(27, dtype=int).reshape(3, 3, 3),
    np.arange(13, dtype=int),
    np.arange(24, dtype=int).reshape(1, 2, 3, 4, 1, 1)]

float_examples = [
    np.arange(27, dtype=float).reshape(3, 3, 3),
    np.arange(13, dtype=float),
    np.arange(24, dtype=float).reshape(1,2,3,4,1,1)]

complex_examples = [
    np.arange(27).reshape(3, 3, 3) * (1+2j),
    np.arange(13) * (1+2j),
    np.arange(24).reshape(1,2,3,4,1,1) * (1+2j)]

string_examples = [
    np.array(["a", "b", "c"])]

examples = int_examples + float_examples + complex_examples + string_examples

@pytest.mark.parametrize("mat", examples)
def test_numpy_serialisation_forwards_backwards(mat: np.ndarray):
    """ End-to-end check for serialisation of numpy arrays"""
    deserialised = numpy_deserialise(numpy_serialise(mat))

    assert mat.dtype == deserialised.dtype
    assert mat.shape == deserialised.shape
    assert np.all(mat.reshape(-1) == deserialised.reshape(-1))


@pytest.mark.parametrize("mat", examples)
def test_numpy_serialisation_bad_size_error(mat):
    """ Numpy serialisations, check that incorrect sizes throw errors """
    serialised = numpy_serialise(mat)
    serialised["data"].append(serialised["data"][-1])

    with pytest.raises(SPWSerialisationError):
        numpy_deserialise(serialised)

@pytest.mark.parametrize("mat", examples)
def test_numpy_serialisation_bad_key_error(mat):
    """ Numpy serialisation, missing key should throw error"""
    serialised = numpy_serialise(mat)

    del serialised["shape"]

    with pytest.raises(SPWSerialisationError):
        numpy_deserialise(serialised)

@pytest.mark.parametrize("mat", float_examples + string_examples)
def test_numpy_serialisation_require_ints(mat):
    """ Numpy serialisation, dtype=int should throw an error with floats or string"""
    serialised = numpy_serialise(mat)

    serialised["dtype"] = np.dtype(int).str

    with pytest.raises(SPWSerialisationError):
        numpy_deserialise(serialised)

