""" Tests for things in util.py """

import pytest

import numpy as np

from util import triple_product_matrix


some_vectors = [np.array(v) for v in [
    [0,0,0],
    [1,0,0],
    [0,1,0],
    [0,0,1],
    [3,4,5]
]]


@pytest.mark.parametrize("a", some_vectors)
@pytest.mark.parametrize("b", some_vectors)
@pytest.mark.parametrize("c", some_vectors)
def test_triple_product_matrix_does_triple_product(a,b,c):
    m = triple_product_matrix(a)
    with_matrix = b @ m @ c
    with_cross = np.dot(a, np.cross(b, c))

    assert with_matrix == pytest.approx(with_cross, abs=1e-9)
