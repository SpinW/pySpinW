import numpy as np
import pytest
from pyspinw import *
from pyspinw.basis import angle_axis_rotation_matrix
from pyspinw.symmetry.operations import SpaceOperation

s2 = np.sqrt(0.5)

expected_off_boundary_points = [
    (0.25, 0.5, 0.5),
    (0.5, 0.75, 0.5),
    (0.25, 0.75, 0.5),
    (0.25, 0.5, 0.5),
]
expected_on_boundary_points = [
    (0, 0.5, 0.5),
    (0.5, 0.5, 0.5),
    (0.5, 0, 0.5),
    (0, 0.5, 0.5),
]
expected_in_spins = [
    (1, 0, 0),
    (0, 1, 0),
    (-s2, -s2, 0),
    (1, 0, 0),
]
expected_round_spins = [
    (0, 1, 0),
    (-s2, -s2, 0),
    (1, 0, 0),
    (0, 1, 0)
]
expected_up_spins = [(0,0,1)]*5

@pytest.mark.parametrize("n_rotations", range(4))
@pytest.mark.parametrize("unit_cell", [
    UnitCell(1,1,1, gamma=120),
    UnitCell(2,2,6, gamma=120)])
@pytest.mark.parametrize("position_or_spins", [True, False])
def test_p6mm_site_operation(unit_cell: UnitCell, n_rotations: int, position_or_spins: bool):
    """ Test that 3-fold rotation about the centre of the cell applies correctly

    The direction of spins on the boundary is complicated,
    because some directions of spin are inconsistent with C4 symmetry
    """

    point_in_off_boundary = LatticeSite(*expected_off_boundary_points[0], 1, 0, 0)
    point_round_off_boundary = LatticeSite(*expected_off_boundary_points[0], 0, 1, 0)
    point_up_off_boundary = LatticeSite(*expected_off_boundary_points[0], 0, 0, 1)

    point_in_on_boundary = LatticeSite(*expected_on_boundary_points[0], 1, 0, 0)
    point_round_on_boundary = LatticeSite(*expected_on_boundary_points[0], 0, 1, 0)
    point_up_on_boundary = LatticeSite(*expected_on_boundary_points[0], 0, 0, 1)

    # Example p6 operation
    rotation_matrix = np.array([
        [ 0, -1, 0],
        [ 1, -1, 0],
        [ 0,  0, 1]
    ])
    translation = np.zeros((3, ))
    operation = SpaceOperation(rotation_matrix, translation)

    for i in range(n_rotations):
        point_in_off_boundary = point_in_off_boundary.symmetry_transformed(operation, unit_cell)
        point_round_off_boundary = point_round_off_boundary.symmetry_transformed(operation, unit_cell)
        point_up_off_boundary = point_up_off_boundary.symmetry_transformed(operation, unit_cell)

        point_in_on_boundary = point_in_on_boundary.symmetry_transformed(operation, unit_cell)
        point_round_on_boundary = point_round_on_boundary.symmetry_transformed(operation, unit_cell)
        point_up_on_boundary = point_up_on_boundary.symmetry_transformed(operation, unit_cell)


    if position_or_spins:
        # Check spatial movement
        assert np.allclose(point_in_off_boundary.ijk, np.array(expected_off_boundary_points[n_rotations]))
        assert np.allclose(point_round_off_boundary.ijk, np.array(expected_off_boundary_points[n_rotations]))
        assert np.allclose(point_up_off_boundary.ijk, np.array(expected_off_boundary_points[n_rotations]))

        assert np.allclose(point_in_on_boundary.ijk, np.array(expected_on_boundary_points[n_rotations]))
        assert np.allclose(point_round_on_boundary.ijk, np.array(expected_on_boundary_points[n_rotations]))
        assert np.allclose(point_up_on_boundary.ijk, np.array(expected_on_boundary_points[n_rotations]))

    else:
        # Check spin directions
        assert np.allclose(point_in_on_boundary.spin_data, np.array(expected_in_spins[n_rotations]))
        assert np.allclose(point_in_off_boundary.spin_data, np.array(expected_in_spins[n_rotations]))

        assert np.allclose(point_round_on_boundary.spin_data, np.array(expected_round_spins[n_rotations]))
        assert np.allclose(point_round_off_boundary.spin_data, np.array(expected_round_spins[n_rotations]))

        assert np.allclose(point_up_on_boundary.spin_data, np.array(expected_up_spins[n_rotations]))
        assert np.allclose(point_up_off_boundary.spin_data, np.array(expected_up_spins[n_rotations]))