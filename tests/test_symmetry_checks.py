import pytest

import numpy as np

from pyspinw.symmetry.symmetry_checking import _symantisym_matrix, _transpose_matrix

def vectorise(mat):
    """ Turn a matrix into vector """
    return mat.reshape(-1)

def unvectorise(vec):
    """ Turn a vector into a matrix """
    return vec.reshape(3,3)


expected = [
    [[1, 0, 0],
     [0, 0, 0],
     [0, 0, 0]],

    [[0, 1, 0],
     [1, 0, 0],
     [0, 0, 0]],

    [[0, 0, 1],
     [0, 0, 0],
     [1, 0, 0]],

    [[0, 0, 0],
     [0, 1, 0],
     [0, 0, 0]],

    [[0, 0, 0],
     [0, 0, 1],
     [0, 1, 0]],

    [[0, 0, 0],
     [0, 0, 0],
     [0, 0, 1]],

    [[0, 0, 0],
     [0, 0, 1],
     [0,-1, 0]],

    [[0, 0,-1],
     [0, 0, 0],
     [1, 0, 0]],

    [[ 0, 1, 0],
     [-1, 0, 0],
     [ 0, 0, 0]],

]

@pytest.mark.parametrize("index_expected_pair", enumerate(expected))
def test_vectorisation_sym_asym(index_expected_pair):
    """ Check that _symantisym makes a matrix of the form

        [[ a, b, c]     [[  0  z -y ]
         [ b, d, e]   +  [ -z  0  x ]
         [ c, e, f]]     [  y -x  0 ]]

    """

    index, exp = index_expected_pair

    v_test = np.zeros((9,))
    v_test[index] = 1

    actual = unvectorise(_symantisym_matrix @ v_test)

    assert np.all(exp == actual)

@pytest.mark.parametrize("index", range(9))
def test_vectorisation_transpose(index):
    v = np.zeros((9,))
    v[index] = 1

    assert np.all(unvectorise(v @ _transpose_matrix) == unvectorise(v).T)
