import numpy as np
import pytest

from pyspinw.site import LatticeSite
from pyspinw.couplinggroup import CouplingGroup, DirectionalityFilter, InDirectionFilter, InPlaneFilter
from pyspinw.coupling import HeisenbergCoupling
from pyspinw.symmetry.unitcell import UnitCell


@pytest.mark.parametrize("lower, upper", [(0, 1), (1, 2), (2, 3), (0, 3)])
def test_simple_coupling_group(lower, upper):
    sites = [LatticeSite(name="X", i=0.5, j=0.5, k=0.5)]

    cell = UnitCell(1, 1, 1)

    group = CouplingGroup(
        name = "test_group",
        min_distance = lower,
        max_distance = upper,
        max_order = None,
        naming_pattern = None,
        coupling_type = HeisenbergCoupling,
        coupling_parameters = {"j": 1.23},
        direction_filter = None)

    couplings = group.couplings(sites, cell)


    # # Can be useful to print this info
    # print()
    # for coupling in couplings:
    #     print(coupling)
    #     print(coupling.distance(cell))

    # Check distances
    for coupling in couplings:
        assert coupling.distance(cell) <= upper, "Coupling distance should not be more than the upper bound"
        assert coupling.distance(cell) >= lower, "Coupling distance should not be less than the lower bound"

        assert coupling.j == 1.23, "J should be set on the couplings"


def test_coupling_group_direction_filtered():
    sites = [LatticeSite(name="X", i=0.5, j=0.5, k=0.5)]

    cell = UnitCell(1, 1, 1)

    group = CouplingGroup(
        name = "test_group",
        min_distance = 0,
        max_distance = 3,
        max_order = None,
        naming_pattern = None,
        coupling_type = HeisenbergCoupling,
        coupling_parameters = {"j": 1.23},
        direction_filter = InDirectionFilter([0,0,1]))

    couplings = group.couplings(sites, cell)

    # Check distances
    for coupling in couplings:
        assert np.sum(np.cross(coupling.vector(cell), [0,0,1])**2) < 1e-10, "Coupling vector should (almost) zero away from z"

    assert len(couplings) > 1, "There should be a few couplings of the specified form"



def test_coupling_group_plane_filtered():
    sites = [LatticeSite(name="X", i=0.5, j=0.5, k=0.5)]

    cell = UnitCell(1, 1, 1)

    group = CouplingGroup(
        name = "test_group",
        min_distance = 0,
        max_distance = 3,
        max_order = None,
        naming_pattern = None,
        coupling_type = HeisenbergCoupling,
        coupling_parameters = {"j": 1.23},
        direction_filter = InPlaneFilter([0,0,1]))

    couplings = group.couplings(sites, cell)

    # Check distances
    for coupling in couplings:
        assert np.sum(np.dot(coupling.vector(cell), [0,0,1])**2) < 1e-10, "Coupling vector should (almost) zero z"

    assert len(couplings) > 1, "There should be a few couplings of the specified form"
