""" Checks for the consistency of ExchangeMatrixConstraints with the checks for specific numerical matrices"""

import pytest

import numpy as np

from pyspinw.exchange import Exchange, DMExchange
from pyspinw import Structure, LatticeSite
from symmetry.symmetry_info_test_cases import cases

test_list = [(structure, site_1, site_2)
                for structure in cases
                for site_1, site_2 in structure.site_pairs()]

rng = np.random.default_rng(667)
@pytest.mark.parametrize("structure, site_1, site_2", test_list)
def test_symmetry_info_and_checks_consistent_positive_results(structure: Structure, site_1: LatticeSite, site_2: LatticeSite):
    """ Generate random values that satisfy the constraints, and check that the check method thinks they're good"""

    exchange_constraints = structure.exchange_constraints(site_1, site_2, do_print=False)

    matrix_values, constraints = exchange_constraints._matrix_form_strings()
    matrix_values = [s.strip() for s in matrix_values]

    for i in range(5):
        # Try five random sets

        lookup = {}
        value_list = []
        for base_string in matrix_values:
            # Zeros should become zeros
            if base_string == '0':
                value_list.append(0.0)
                continue

            # remove any leading - for testing
            if base_string.startswith("-"):
                string = base_string[1:]
            else:
                string = base_string

            if string in lookup:
                value = lookup[string]
            else:
                value = 2*rng.random() - 1
                lookup[string] = value

            if base_string.startswith("-"):
                value *= -1

            value_list.append(value)

        # Create exchanges matching the values
        l = value_list
        symmetric_exchange_matrix = np.array([
            [l[0], l[1], l[2]],
            [l[1], l[3], l[4]],
            [l[2], l[4], l[5]]])
        symmetric_exchange = Exchange(site_1, site_2, exchange_matrix=symmetric_exchange_matrix, cell_offset=(0,0,0))

        # Check that they are considered correct

        assert symmetric_exchange.obeys_symmetry(structure), (f"Symmetric exchange should pass check,"
                                                              f" matrix is {symmetric_exchange_matrix}")


        # Only check DM if not zero
        if not (l[6] == 0 and l[7]==0 and l[8]==0):
            dm_exchange = DMExchange(site_1, site_2, d_x=l[6],d_y=l[7], d_z=l[8], cell_offset=(0,0,0))
            assert dm_exchange.obeys_symmetry(structure), f"DM Exchange should pass check, vector is {l[6:]}"


            # Check the combined exchange matrix
            full_exchange = Exchange(site_1, site_2,
                                     exchange_matrix=symmetric_exchange_matrix + dm_exchange.exchange_matrix,
                                     cell_offset=(0,0,0))

            assert full_exchange.obeys_symmetry(structure), ("Full asymmetric exchange should pass check, matrix is "
                                                             f"{full_exchange.exchange_matrix}")


