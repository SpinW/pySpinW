""" Tests for the various symmetry checking routines, this should give an
opportunity to check the related xyz/rlu transforms """

from pyspinw import *
import pytest

from symmetry.structure_test_cases import cases

@pytest.mark.parametrize("case", cases)
@pytest.mark.parametrize("neighbour", range(1, 4))
def test_heisenberg_always_allowed_exchange_obeys_symmetry(case, neighbour):
    """ Heisenberg exchanges should always be allowed, check `Exchange.obeys_symmetry` does this"""

    structure = case()

    for site_1 in structure.sites:
        for site_2, offset in structure.neighbours(site_1, n=neighbour):
            ex = HeisenbergExchange(site_1, site_2, cell_offset=offset, j=-1)
            assert ex.obeys_symmetry(structure)


@pytest.mark.parametrize("case", cases)
@pytest.mark.parametrize("neighbour", range(1, 4))
def test_heisenberg_always_allowed_structure_exchange_constraints(case, neighbour):
    """ Heisenberg exchanges should always be allowed, check `Structure.exchange_constraints` does this"""

    structure = case()

    for site_1 in structure.sites:
        for site_2, offset in structure.neighbours(site_1, n=neighbour):
            constraints = structure.exchange_constraints(site_1, site_2)

            assert constraints.text_summary()


