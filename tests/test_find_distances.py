""" Tests for the algorithm to find vectors within a crystal """
import numpy as np
import pytest

from pyspinw.lattice_distances import find_relative_positions
from pyspinw.symmetry.unitcell import UnitCell

simple_unit_cell = UnitCell(1,1,1)
squashed_unit_cell = UnitCell(1, 10, 10)

@pytest.mark.parametrize("distance", [0.0, 0.1, 0.4])
def test_single_distance_base_cell(distance):
    difference = np.array([-distance, 0, 0])
    output = find_relative_positions(difference, unit_cell_transform=simple_unit_cell._xyz, max_distance=0.5)

    assert len(output.cell_indices.reshape(-1)) == 3, "Should be one 3-vector"
    assert np.all(output.distances == distance), "Should have same distance"
    assert np.all(output.cell_indices == np.array([0,0,0])), "Should have (0,0,0) cell offset"


@pytest.mark.parametrize("distance", [0.6, 0.9])
def test_single_distance_neighbour_cell(distance):
    difference = np.array([-distance, 0, 0])
    output = find_relative_positions(difference, unit_cell_transform=simple_unit_cell._xyz, max_distance=0.5)

    assert len(output.cell_indices.reshape(-1)) == 3, "Should be one 3-vector"
    assert np.all(output.distances == 1-distance), "Should have 1 - distance"
    assert np.all(output.cell_indices == np.array([1,0,0])), "Should have (1,0,0) cell offset"

@pytest.mark.parametrize("distance", [0.1, 0.4, 0.6, 0.9])
def test_double_distance(distance):
    difference = np.array([-distance, 0, 0])
    output = find_relative_positions(difference, unit_cell_transform=squashed_unit_cell._xyz, max_distance=1.0)

    assert len(output.cell_indices.reshape(-1)) == 6, "Should be two 3-vectors"
    assert np.all(np.logical_or(output.distances == distance, output.distances == 1-distance)), "Distances should be either distance or 1-distance"

