""" Tests for supercells """
import numpy as np
from pyspinw.cell_offsets import CellOffset
from pyspinw.site import LatticeSite
from pyspinw.symmetry.supercell import (PropagationVector, CommensuratePropagationVector, RotationTransform,
                                        TrivialSupercell, TransformationSupercell, SummationSupercell,
                                        RotationSupercell)

def test_minimal_cell_size_1():
    v1 = CommensuratePropagationVector(i=1/2, j=1, k=1)
    v2 = CommensuratePropagationVector(i=1/3, j=1, k=1)

    supercell = SummationSupercell(propagation_vectors=[v1, v2])

    i,j,k = supercell.cell_size()

    assert i == 6
    assert j == 1
    assert k == 1

    assert supercell.n_cells() == 6
    assert supercell.scaling == (6, 1, 1)

def test_antiferromagnet_chain():
    k = CommensuratePropagationVector(i=1/2, j=1, k=1)
    x = LatticeSite(0, 0, 0, 0, 1, 0, name="X")
    delta = CellOffset(1, 0, 0)

    supercell = SummationSupercell(propagation_vectors=[k])
    
    moment = supercell.moment(x, delta)

    assert np.allclose(moment, np.array([0, -1, 0]))

def test_triangular_antiferromagnet():
    k = CommensuratePropagationVector(i=1/3, j=1, k=1)
    x = LatticeSite(0, 0, 0, 0, 1, 0, name="X")
    delta = CellOffset(1, 0, 0)
    rotmat = RotationTransform(axis=[0, 0, 1])

    supercell = TransformationSupercell([(k, rotmat)])
    moment = supercell.moment(x, delta)
    assert np.allclose(moment, np.array([-np.sqrt(3)/2, -0.5, 0]))

    # Now try it with the "incommensurate" supercell
    supercell = RotationSupercell(perpendicular=[0, 0, 1], propagation_vector=[1/3, 0, 0])
    moment = supercell.moment(x, delta)
    assert np.allclose(moment, np.array([-np.sqrt(3)/2, -0.5, 0]))

def test_rotationsupercell_approximant():
    supercell = RotationSupercell(perpendicular=[0, 0, 1], propagation_vector=[1/3, 0, 0])
    approx = supercell.approximant()

    assert isinstance(approx, TransformationSupercell)
    assert len(approx._transforms) == 1
    assert len(approx._propagation_vectors) == 1
    assert np.allclose(approx._propagation_vectors[0]._vector[0], 1/3)
    assert isinstance(approx._transforms[0][1], RotationTransform)
    assert np.allclose(approx._transforms[0][1]._axis, [0, 0, 1])
