import pytest

from pyspinw.exchange import HeisenbergExchange
from pyspinw.exchangegroup import DirectionalityFilter, InDirectionFilter, BiDirectionFilter
from pyspinw.interface import generate_exchanges
from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell

directions = [[1,0,0],
              [0,1,0],
              [0,0,1]]

@pytest.mark.parametrize("direction", directions)
def test_heisenberg_ferromagnet(direction):
    """ Test that a simple Heisenberg ferromagnet only gives one exchange"""
    sites = [LatticeSite(0,0,0,0,0,1, name="X")]

    exchanges = generate_exchanges(sites, UnitCell(1, 1, 1), max_distance=1,
                                   direction_filter=BiDirectionFilter(direction), j=-1)

    assert len(exchanges) == 1


def test_heisenberg_antiferromagnet_supercell():

    sites = [LatticeSite(0,0,0,0,0,1, name="X"),
             LatticeSite(0.5,0,0,0,0,-1, name="Y")]

    exchanges = generate_exchanges(sites, UnitCell(1, 1, 1), max_distance=0.6,
                                   direction_filter=BiDirectionFilter([1,0,0]), j=1)

    assert len(exchanges) == 2

def test_heisenberg_antiferromagnet_single_cell():
    sites = [LatticeSite(0, 0, 0, 0, 0, 1, name="X")]

    exchanges = generate_exchanges(sites, UnitCell(1, 1, 1), max_distance=1.1,
                                   direction_filter=BiDirectionFilter([1, 0, 0]), j=1)

    assert len(exchanges) == 1


def test_hexagonal_cell():
    sites = [LatticeSite(0, 0, 0, 0, 0, 1, name="X")]

    exchanges = generate_exchanges(sites, UnitCell(1, 1, 2, gamma=120), max_distance=1.1, j=1)

    assert len(exchanges) == 3


def test_cubic_cell():
    sites = [LatticeSite(0, 0, 0, 0, 0, 1, name="X")]

    # Nearest neighbours (1 in each x, y, z directions)
    exchanges = generate_exchanges(sites, UnitCell(1, 1, 1), bond=1, j=1)
    assert len(exchanges) == 3

    # Next-nearest neighbours (1 each in +x+y, +x-y, +x+z, +x-z, +y+z, +y-z directions)
    exchanges = generate_exchanges(sites, UnitCell(1, 1, 1), bond=2, j=1)
    assert len(exchanges) == 6

    # 3rd-nearest neighbours (1 in each x+y+z, x-y+z, x+y-z, x-y-z directions)
    exchanges = generate_exchanges(sites, UnitCell(1, 1, 1), bond=3, j=1)
    assert len(exchanges) == 4

