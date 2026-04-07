import numpy as np
import pytest

from pyspinw.site import LatticeSite
from pyspinw.exchangegroup import ExchangeGroup, DirectionalityFilter, InDirectionFilter, InPlaneFilter
from pyspinw.exchange import HeisenbergExchange
from pyspinw.symmetry.unitcell import UnitCell


@pytest.mark.parametrize("lower, upper", [(0, 1), (1, 2), (2, 3), (0, 3)])
def test_simple_coupling_group(lower, upper):
    sites = [LatticeSite(name="X", i=0.5, j=0.5, k=0.5)]

    cell = UnitCell(1, 1, 1)

    group = ExchangeGroup(
        name = "test_group",
        bond = 0,
        min_distance = lower,
        max_distance = upper,
        max_order = None,
        naming_pattern = None,
        exchange_type= HeisenbergExchange,
        coupling_parameters = {"j": 1.23},
        direction_filter = None)

    couplings = group.exchanges(sites, cell)


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

    group = ExchangeGroup(
        name = "test_group",
        bond = 0,
        min_distance = 0,
        max_distance = 3,
        max_order = None,
        naming_pattern = None,
        exchange_type= HeisenbergExchange,
        coupling_parameters = {"j": 1.23},
        direction_filter = InDirectionFilter([0,0,1]))

    couplings = group.exchanges(sites, cell)

    # Check distances
    for coupling in couplings:
        assert np.sum(np.cross(coupling.vector(cell), [0,0,1])**2) < 1e-10, "Coupling vector should (almost) zero away from z"

    assert len(couplings) > 1, "There should be a few couplings of the specified form"



def test_coupling_group_plane_filtered():
    sites = [LatticeSite(name="X", i=0.5, j=0.5, k=0.5)]

    cell = UnitCell(1, 1, 1)

    group = ExchangeGroup(
        name = "test_group",
        bond = 0,
        min_distance = 0,
        max_distance = 3,
        max_order = None,
        naming_pattern = None,
        exchange_type= HeisenbergExchange,
        coupling_parameters = {"j": 1.23},
        direction_filter = InPlaneFilter([0,0,1]))

    couplings = group.exchanges(sites, cell)

    # Check distances
    for coupling in couplings:
        assert np.sum(np.dot(coupling.vector(cell), [0,0,1])**2) < 1e-10, "Coupling vector should (almost) zero z"

    assert len(couplings) > 1, "There should be a few couplings of the specified form"


def test_coupling_bond_index():
    sites = [LatticeSite(name="X", i=0., j=0., k=0.)]

    cell = UnitCell(1, 1, 2, gamma=120.)

    group = ExchangeGroup(
        name = "test_group",
        bond = 1,
        min_distance = 0,
        max_distance = 0,
        max_order = None,
        naming_pattern = None,
        exchange_type= HeisenbergExchange,
        coupling_parameters = {"j": 1.23},
        direction_filter = None)

    couplings = group.exchanges(sites, cell)

    # Check distances
    for coupling in couplings:
        assert np.allclose(coupling.distance(cell), 1.0), "Coupling distance should (almost) unity"

    # Only 3 couplings because we only have the forward directions, [100], [110] and [010]
    assert len(couplings) == 3, "There should be 3 couplings of the specified form"
