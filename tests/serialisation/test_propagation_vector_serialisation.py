import pytest

from pyspinw.symmetry.supercell import PropagationVector, CommensuratePropagationVector

some_numbers = [0.1, 1/3, 4/7, 18/21, 0]

@pytest.mark.parametrize("i", some_numbers)
@pytest.mark.parametrize("j", some_numbers)
@pytest.mark.parametrize("k", some_numbers)
def test_incommensurate(i,j,k):
    vector = PropagationVector(i,j,k)
    serialised = vector.serialise()
    deserialised = PropagationVector.deserialise(serialised)

    assert isinstance(deserialised, PropagationVector)
    assert not isinstance(deserialised, CommensuratePropagationVector)

    assert deserialised.i == vector.i
    assert deserialised.j == vector.j
    assert deserialised.k == vector.k


@pytest.mark.parametrize("i", some_numbers)
@pytest.mark.parametrize("j", some_numbers)
@pytest.mark.parametrize("k", some_numbers)
def test_commensurate(i, j, k):
    vector = CommensuratePropagationVector(i, j, k)
    serialised = vector.serialise()
    deserialised = PropagationVector.deserialise(serialised)

    assert isinstance(deserialised, CommensuratePropagationVector)

    assert deserialised.i == vector.i
    assert deserialised.j == vector.j
    assert deserialised.k == vector.k
