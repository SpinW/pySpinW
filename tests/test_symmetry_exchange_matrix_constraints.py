""" Tests for ExchangeMatrixConstraints """

import pytest

from pyspinw import *

from pyspinw.symmetry.symmetry_checking import ExchangeMatrixConstraints

from symmetry.exchange_and_anisotropy_matrix_constraints_object_cases import ops, sym, asym, string_to_ops

cell = UnitCell(1,1,1)

@pytest.mark.parametrize("operations_string, symmetric, asymmetric", zip(ops, sym, asym))
def test_exchange_matrix_constraints_object_gives_correct_results(operations_string, symmetric, asymmetric):
    """ Check the output strings match the expected output from Mathematica"""
    operations = string_to_ops(operations_string)
    transforms = [op.point_operation_in_cartesian(cell) for op in operations]
    constraints = ExchangeMatrixConstraints(transforms)

    expressions, _ = constraints._matrix_form_strings()
    expressions = [expr.strip() for expr in expressions]

    # What should the list of strings look like, we can assume that we don't need to worry about constraints,
    #  or complicated expressions, because they don't appear in any of actual test cases (this could change!)

    test_values = [entry
                    for i, entries in enumerate(symmetric)
                    for j, entry in enumerate(entries)
                        if i <= j] + asymmetric

    print(expressions)
    print(test_values)

    lookup = {}
    for actual, expected in zip(expressions, test_values):
        if expected == 0:
            assert actual == '0', f"Zero entries should be the same, (got '{actual}' not '0')"

        else:
            assert actual != '0', f"Expected entry matching {expected} to not be zero"
            if expected in lookup:
                assert lookup[expected] == actual, (f"Same expected number should mean same letter "
                                                    f"(lookup[{expected}] = {lookup[expected]}, not {actual})")
            else:
                lookup[expected] = actual

    # check for negations
    for key in lookup:
        if key > 0:
            if -key in lookup:
                assert lookup[key] == "-"+lookup[-key] or lookup[-key] == "-" + lookup[key], ("Negated key should"
                                                                                              "match negated entry")


