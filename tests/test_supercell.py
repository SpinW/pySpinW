""" Tests for supercells """
from pyspinw.supercell import CommensuratePropagationVector, CommensurateSupercell


def test_minimal_cell_size_1():
    v1 = CommensuratePropagationVector(i=1/2, j=1, k=1)
    v2 = CommensuratePropagationVector(i=1/3, j=1, k=1)

    supercell = CommensurateSupercell(propagation_vectors=[v1, v2])

    i,j,k = supercell.minimal_supercell()

    assert i == 6
    assert j == 1
    assert k == 1

