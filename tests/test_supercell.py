""" Tests for supercells """
import numpy as np
import pytest
from pyspinw.cell_offsets import CellOffset
from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structures import Structure
from pyspinw.symmetry.supercell import TransformationSupercell, CommensuratePropagationVector, RotationTransform, \
    TiledSupercell, SummationSupercell, RotationSupercell


def test_trivial_supercell():
    """ Test that trivial supercell works as it should """

    unit_cell = UnitCell(1, 1, 1)

    x = LatticeSite(0, 0, 0, 1, 0, 0, name="X")

    supercell = TiledSupercell((1, 2, 3))

    structure = Structure([x], unit_cell=unit_cell, spacegroup=None, supercell=supercell)

    expanded = structure.expand()

    assert len(expanded.sites) == 1 * 2 * 3

    assert supercell.n_components() == 1


@pytest.mark.parametrize("n", [3, 4, 5])
def test_transform_supercell(n):
    """ Check that transformation supercell works as it should """
    unit_cell = UnitCell(1,1,1)

    x = LatticeSite(0,0,0,1,0,0, name="X")

    propagation_vector = CommensuratePropagationVector(0,0,1/n)

    supercell = TransformationSupercell([(propagation_vector, RotationTransform([0,0,1]))])

    structure = Structure([x], unit_cell=unit_cell, spacegroup=None, supercell=supercell)

    expanded = structure.expand()

    for site1, site2 in zip(expanded.sites, expanded.sites[1:]):
        assert np.isclose(np.dot(site1.base_spin, site2.base_spin), np.cos(2 * np.pi / n))

    assert supercell.n_components() == 1

@pytest.mark.parametrize("n", [3, 4, 5])
def test_summation_supercell(n):
    """ Check that transformation supercell works as it should """
    unit_cell = UnitCell(1,1,1)

    x = LatticeSite(0, 0, 0, supercell_spins=[[1, 0, 0], [0, 1, 0]], name="X")

    propagation_vector_1 = CommensuratePropagationVector(0,0, 1/n) # Cosine
    propagation_vector_2 = CommensuratePropagationVector(0,0, 1/n, phase=np.pi/2) # -Sine

    supercell = SummationSupercell([propagation_vector_1, propagation_vector_2])

    structure = Structure([x], unit_cell=unit_cell, spacegroup=None, supercell=supercell)

    expanded = structure.expand()

    for site1, site2 in zip(expanded.sites, expanded.sites[1:]):
        assert np.isclose(np.dot(site1.base_spin, site2.base_spin), np.cos(2 * np.pi / n))

    assert supercell.n_components() == 2

def test_antiferromagnet_chain():
    k = CommensuratePropagationVector(i=1/2, j=1, k=1)
    x = LatticeSite(0, 0, 0, 0, 1, 0, name="X")
    delta = CellOffset(1, 0, 0)

    supercell = SummationSupercell(propagation_vectors=[k])
    
    moment = supercell.spin(x, delta)

    assert np.allclose(moment, np.array([0, -1, 0]))

def test_triangular_antiferromagnet():
    k = CommensuratePropagationVector(i=1/3, j=1, k=1)
    x = LatticeSite(0, 0, 0, 0, 1, 0, name="X")
    delta = CellOffset(1, 0, 0)
    rotmat = RotationTransform(axis=[0, 0, 1])

    supercell = TransformationSupercell([(k, rotmat)])
    moment = supercell.spin(x, delta)
    assert np.allclose(moment, np.array([-np.sqrt(3)/2, -0.5, 0]))

    # Now try it with the "incommensurate" supercell
    supercell = RotationSupercell(perpendicular=[0, 0, 1], propagation_vector=[1/3, 0, 0])
    moment = supercell.spin(x, delta)
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
