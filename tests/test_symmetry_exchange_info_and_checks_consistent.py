""" Checks for the consistency of AnisotropyMatrixConstraints with the checks for specific numerical matrices"""
from collections import defaultdict

import pytest

import numpy as np

from pyspinw import Structure, LatticeSite, Anisotropy
from symmetry.symmetry_info_test_cases import cases

test_list = [(structure, site)
                for structure in cases
                for site in structure.sites]

rng = np.random.default_rng(667)

@pytest.mark.parametrize("structure, site", test_list, ids=[f"{x[1].name}" for x in test_list])
def test_anisotropy_info_and_checks_consistent_positive_results(
        structure: Structure, site: LatticeSite):
    """ Generate random values that satisfy the constraints, and check that the check method thinks they're good"""

    anisotropy_constraints = structure.anisotropy_constraints(site, do_print=False)

    matrix_values, constraints = anisotropy_constraints._matrix_form_strings()
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

        # Create anisotropies matching the values
        l = value_list
        symmetric_anisotropy_matrix = np.array([
            [l[0], l[1], l[2]],
            [l[1], l[3], l[4]],
            [l[2], l[4], l[5]]])
        anisotropy = Anisotropy(site, anisotropy_matrix=symmetric_anisotropy_matrix)

        # Check that they are considered correct

        assert anisotropy.obeys_symmetry(structure), (f"Anisotropy should pass check,"
                                                              f" matrix is {symmetric_anisotropy_matrix}")


@pytest.mark.parametrize("structure, site", test_list, ids=[f"{x[1].name}" for x in test_list])
def test_anisotropy_info_and_checks_consistent_negative_results_zeros(
        structure: Structure, site: LatticeSite):
    """ Generate random values that don't satisfy the constraints by having non-zero values where there should be zeros

    Check that the check is false!
    """

    anisotropy_constraints = structure.anisotropy_constraints(site, do_print=False)

    matrix_values, constraints = anisotropy_constraints._matrix_form_strings()
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

        # Create anisotropy matching the values
        l = value_list
        matrix = np.array([
            [l[0], l[1], l[2]],
            [l[1], l[3], l[4]],
            [l[2], l[4], l[5]]])
        anisotropy = Anisotropy(site, anisotropy_matrix=matrix)

        # Check that they are considered wrong

        assert not anisotropy.obeys_symmetry(structure), (f"Anisotropy should fail check,"
                                                              f" matrix is {matrix}")




@pytest.mark.parametrize("structure, site", test_list, ids=[f"{x[1].name}" for x in test_list])
def test_anisotropy_info_and_checks_consistent_negative_results_break_equality(
        structure: Structure, site: LatticeSite):
    """ Generate random values that don't satisfy the constraints by having different values where they should be same

    Check that the check is false!
    """

    anisotropy_constraints = structure.anisotropy_constraints(site, do_print=False)

    matrix_values, constraints = anisotropy_constraints._matrix_form_strings()
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

        # Create anisotropy matching the values
        l = value_list
        matrix = np.array([
            [l[0], l[1], l[2]],
            [l[1], l[3], l[4]],
            [l[2], l[4], l[5]]])
        anisotropy = Anisotropy(site, anisotropy_matrix=matrix)

        # Check that they are considered wrong

        assert not anisotropy.obeys_symmetry(structure), (f"Anisotropy should fail check,"
                                                              f" matrix is {matrix}")
