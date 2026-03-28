""" Tests for supercells """
import numpy as np
import pytest

from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structures import Structure
from pyspinw.symmetry.supercell import TransformationSupercell, CommensuratePropagationVector, RotationTransform, \
    TrivialSupercell, SummationSupercell


def test_trivial_supercell():
    """ Test that trivial supercell works as it should """

    unit_cell = UnitCell(1, 1, 1)

    x = LatticeSite(0, 0, 0, 1, 0, 0, name="X")

    supercell = TrivialSupercell((1, 2, 3))

    structure = Structure([x], unit_cell=unit_cell, spacegroup=None, supercell=supercell)

    expanded = structure.expand()

    assert len(expanded.sites) == 1*2*3

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
        assert np.isclose(np.dot(site1.base_moment, site2.base_moment), np.cos(2*np.pi/n))

    assert supercell.n_components() == 1

@pytest.mark.parametrize("n", [3, 4, 5])
def test_summation_supercell(n):
    """ Check that transformation supercell works as it should """
    unit_cell = UnitCell(1,1,1)

    x = LatticeSite(0,0,0,supercell_moments=[[1,0,0],[0,1,0]], name="X")

    propagation_vector_1 = CommensuratePropagationVector(0,0, 1/n) # Cosine
    propagation_vector_2 = CommensuratePropagationVector(0,0, 1/n, phase=np.pi/2) # -Sine

    supercell = SummationSupercell([propagation_vector_1, propagation_vector_2])

    structure = Structure([x], unit_cell=unit_cell, spacegroup=None, supercell=supercell)

    expanded = structure.expand()

    for site1, site2 in zip(expanded.sites, expanded.sites[1:]):
        assert np.isclose(np.dot(site1.base_moment, site2.base_moment), np.cos(2*np.pi/n))

    assert supercell.n_components() == 2

