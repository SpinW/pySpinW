""" Tests for things in basis.py """

import pytest
import numpy as np

from pyspinw.basis import find_aligned_basis

#
# find_aligned_basis tests
#

# Random, but deterministic, vectors for testing
input_lengths = [0, 1, 3, 100] # Definitely need cases of 0, 1 and 3 for potential issues around array shapes
rng = np.random.default_rng(1984)
test_vectors = [rng.normal(size=(n, 3)) for n in input_lengths]

test_vectors.append(np.eye(3)) # Basis aligned unit vectors
test_vectors.append(10*np.eye(3)) # Basis aligned non-unit vectors

acceptable_tolerance = 10*np.finfo(float).eps

zeros = [np.zeros((10, 3), dtype=float)] # Sometimes we want to check length zero input, sometimes not

@pytest.mark.parametrize("vectors", test_vectors + zeros)
def test_find_aligned_basis_normality(vectors):
    """ Check length of all basis vectors is 1 within tolerance"""
    for basis_vector in find_aligned_basis(vectors):
        lengths = np.sqrt(np.sum(basis_vector**2, axis=1))
        errors = lengths - 1
        assert np.all(np.abs(errors) < acceptable_tolerance)


@pytest.mark.parametrize("vectors", test_vectors + zeros)
def test_find_aligned_basis_orthogonality(vectors):
    """ Check that the given basis vectors are orthogonal """
    e1, e2, e3 = find_aligned_basis(vectors)

    d12 = np.sum(e1 * e2, axis=1)
    d23 = np.sum(e2 * e3, axis=1)
    d31 = np.sum(e3 * e1, axis=1)

    assert np.all(np.abs(d12) < acceptable_tolerance)
    assert np.all(np.abs(d23) < acceptable_tolerance)
    assert np.all(np.abs(d31) < acceptable_tolerance)


@pytest.mark.parametrize("vectors", test_vectors)
def test_find_aligned_basis_matches(vectors):
    """ Check that the aligned basis vectors are in fact aligned

    use x.y= |x| if y is a unit vector aligned with x, i.e. x.y = |x||y| for parallel vectors

    """
    e1, _, _ = find_aligned_basis(vectors)

    vector_lengths = np.sqrt(np.sum(vectors**2, axis=1))

    dot_product = np.sum(e1 * vectors, axis=1)

    assert np.all(np.abs(dot_product - vector_lengths) < acceptable_tolerance)

def test_find_aligned_basis_zero_is_z():
    """ Check that a sensible error is thrown if any vector """
    e1, e2, e3 = find_aligned_basis(zeros[0])

    assert np.all(e1[:, 0] == 0)
    assert np.all(e1[:, 1] == 0)
    assert np.all(e1[:, 2] == 1)

    assert np.all(e2[:, 0] == 0)
    assert np.all(e2[:, 1] == -1)
    assert np.all(e2[:, 2] == 0)

    assert np.all(e3[:, 0] == 1)
    assert np.all(e3[:, 1] == 0)
    assert np.all(e3[:, 2] == 0)
