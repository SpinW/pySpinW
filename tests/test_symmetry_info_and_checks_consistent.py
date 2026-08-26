""" Checks for the consistency of ExchangeMatrixConstraints with the checks for specific numerical matrices"""
from collections import defaultdict

import pytest

import numpy as np

from pyspinw.exchange import Exchange, DMExchange
from pyspinw import Structure, LatticeSite
from symmetry.symmetry_info_test_cases import cases

test_list = [(structure, site_1, site_2)
                for structure in cases
                for site_1, site_2 in structure.site_pairs()]

rng = np.random.default_rng(667)

@pytest.mark.parametrize("structure, site_1, site_2", test_list,
                         ids=[f"{x[1].name} -> {x[2].name}" for x in test_list])
def test_symmetry_info_and_checks_consistent_positive_results(
        structure: Structure, site_1: LatticeSite, site_2: LatticeSite):
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


@pytest.mark.parametrize("structure, site_1, site_2", test_list,
                         ids=[f"{x[1].name} -> {x[2].name}" for x in test_list])
def test_symmetry_info_and_checks_consistent_negative_results_zeros(
        structure: Structure, site_1: LatticeSite, site_2: LatticeSite):
    """ Generate random values that don't satisfy the constraints by having non-zero values where there should be zeros

    Check that the check is false!
    """

    exchange_constraints = structure.exchange_constraints(site_1, site_2, do_print=False)

    matrix_values, constraints = exchange_constraints._matrix_form_strings()
    matrix_values = [s.strip() for s in matrix_values]

    # We can't break zero cases unless there are zero entries
    if not np.any([value == '0' for value in matrix_values]):
        return

    for i in range(5):
        # Try five random sets

        zeros = np.where([value == '0' for value in matrix_values])[0]
        index = zeros[rng.integers(len(zeros))]

        tampered_matrix_values = matrix_values.copy()
        tampered_matrix_values[int(index)] = "REPLACE" # This should be replaced by a number

        lookup = {}
        value_list = []
        for base_string in tampered_matrix_values:
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

        # Create exchange matching the values
        l = value_list
        x = l[6]
        y = l[7]
        z = l[8]
        exchange_matrix = np.array([
            [l[0], l[1] + z, l[2] - y],
            [l[1]-z, l[3], l[4] + x],
            [l[2]+y, l[4]-x, l[5]]])
        exchange = Exchange(site_1, site_2, exchange_matrix=exchange_matrix, cell_offset=(0,0,0))

        # Check that they are considered wrong

        assert not exchange.obeys_symmetry(structure), (f"Exchange should fail check,"
                                                              f" matrix is {exchange_matrix}")



@pytest.mark.parametrize("structure, site_1, site_2", test_list,
                         ids=[f"{x[1].name} -> {x[2].name}" for x in test_list])
def test_symmetry_info_and_checks_consistent_negative_results_break_equality(
        structure: Structure, site_1: LatticeSite, site_2: LatticeSite):
    """ Generate random values that don't satisfy the constraints by having different values where they should be same

    Check that the check is false!
    """

    exchange_constraints = structure.exchange_constraints(site_1, site_2, do_print=False)

    matrix_values, constraints = exchange_constraints._matrix_form_strings()
    matrix_values = [s.strip() for s in matrix_values]

    # We need to have at least two entries that are the same for this to work
    lookup = defaultdict(list)
    for index, value in enumerate(matrix_values):
        unnegative = value[1:] if value.startswith("-") else value
        lookup[unnegative].append(index)

    if not np.any([len(x) > 1 for x in lookup.values()]):
        return

    # Find potential places to make a replacement
    potential_replacements = [index
                              for indices in lookup.values()
                              for index in indices
                              if len(indices) > 1]


    for i in range(5):
        # Try five random sets


        index = potential_replacements[rng.integers(len(potential_replacements))]

        tampered_matrix_values = matrix_values.copy()
        tampered_matrix_values[index] = "REPLACE" # This should be replaced by a number



        lookup = {}
        value_list = []
        for base_string in tampered_matrix_values:
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

        # Create exchange matching the values
        l = value_list
        x = l[6]
        y = l[7]
        z = l[8]
        exchange_matrix = np.array([
            [l[0], l[1] + z, l[2] - y],
            [l[1]-z, l[3], l[4] + x],
            [l[2]+y, l[4]-x, l[5]]])
        exchange = Exchange(site_1, site_2, exchange_matrix=exchange_matrix, cell_offset=(0,0,0))

        # Check that they are considered wrong

        assert not exchange.obeys_symmetry(structure), (f"Exchange should fail check,"
                                                              f" matrix is {exchange_matrix}")
