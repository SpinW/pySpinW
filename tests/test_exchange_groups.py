import numpy as np
import pytest

from pyspinw.site import LatticeSite
from pyspinw.exchangegroup import ExchangeGroup, DirectionalityFilter, InDirectionFilter, InPlaneFilter
from pyspinw.exchange import HeisenbergExchange
from pyspinw.symmetry.unitcell import UnitCell


@pytest.mark.parametrize("lower, upper", [(0, 1), (1, 2), (2, 3), (0, 3)])
def test_simple_exchange_group(lower, upper):
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
        exchange_parameters = {"j": 1.23},
        direction_filter = None)

    exchanges = group.exchanges(sites, cell)


    # # Can be useful to print this info
    # print()
    # for exchange in exchanges:
    #     print(exchange)
    #     print(exchange.distance(cell))

    # Check distances
    for exchange in exchanges:
        assert exchange.distance(cell) <= upper, "Exchange distance should not be more than the upper bound"
        assert exchange.distance(cell) >= lower, "Exchange distance should not be less than the lower bound"

        assert exchange.j == 1.23, "J should be set on the exchanges"


def test_exchange_group_direction_filtered():
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
        exchange_parameters = {"j": 1.23},
        direction_filter = InDirectionFilter([0,0,1]))

    exchanges = group.exchanges(sites, cell)

    # Check distances
    for exchange in exchanges:
        assert np.sum(np.cross(exchange.vector(cell), [0,0,1])**2) < 1e-10, "Exchange vector should (almost) zero away from z"

    assert len(exchanges) > 1, "There should be a few exchanges of the specified form"



def test_exchange_group_plane_filtered():
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
        exchange_parameters = {"j": 1.23},
        direction_filter = InPlaneFilter([0,0,1]))

    exchanges = group.exchanges(sites, cell)

    # Check distances
    for exchange in exchanges:
        assert np.sum(np.dot(exchange.vector(cell), [0,0,1])**2) < 1e-10, "Exchange vector should (almost) zero z"

    assert len(exchanges) > 1, "There should be a few exchanges of the specified form"


def test_exchange_bond_index():
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
        exchange_parameters = {"j": 1.23},
        direction_filter = None)

    exchanges = group.exchanges(sites, cell)

    # Check distances
    for exchange in exchanges:
        assert np.allclose(exchange.distance(cell), 1.0), "Exchange distance should (almost) unity"

    # Only 3 exchanges because we only have the forward directions, [100], [110] and [010]
    assert len(exchanges) == 3, "There should be 3 exchanges of the specified form"
