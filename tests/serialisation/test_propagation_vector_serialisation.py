import numpy as np
import pytest

from pyspinw.symmetry.supercell import PropagationVector, CommensuratePropagationVector

some_numbers = [0.1, 1/3, 4/7, 18/21, 0]
some_irrational_numbers = [np.pi, np.sqrt(2), np.exp(50), np.exp(-50)]

@pytest.mark.parametrize("i", some_numbers + some_irrational_numbers)
@pytest.mark.parametrize("j", some_numbers + some_irrational_numbers)
@pytest.mark.parametrize("k", some_numbers + some_irrational_numbers)
@pytest.mark.parametrize("phase", [0.0, np.pi/2, 0.1])
def test_incommensurate(i,j,k, phase):
    vector = PropagationVector(i,j,k, phase)
    serialised = vector.serialise()
    deserialised = PropagationVector.deserialise(serialised)

    assert isinstance(deserialised, PropagationVector)
    assert not isinstance(deserialised, CommensuratePropagationVector)

    assert deserialised.i == vector.i
    assert deserialised.j == vector.j
    assert deserialised.k == vector.k
    assert deserialised.phase == vector.phase


@pytest.mark.parametrize("i", some_numbers)
@pytest.mark.parametrize("j", some_numbers)
@pytest.mark.parametrize("k", some_numbers)
@pytest.mark.parametrize("phase", [0.0, np.pi/2, 0.1])
def test_commensurate(i, j, k, phase):
    vector = CommensuratePropagationVector(i, j, k, phase)
    serialised = vector.serialise()
    deserialised = PropagationVector.deserialise(serialised)

    assert isinstance(deserialised, CommensuratePropagationVector)

    assert deserialised.i == vector.i
    assert deserialised.j == vector.j
    assert deserialised.k == vector.k
    assert deserialised.phase == vector.phase
