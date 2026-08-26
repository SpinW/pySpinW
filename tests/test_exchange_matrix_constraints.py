""" Tests for ExchangeMatrixConstraints """

import pytest

from pyspinw import *

from pyspinw.symmetry.symmetry_checking import ExchangeMatrixConstraints

from symmetry.exchange_matrix_contraints_object_cases import ops, sym, asym, string_to_ops

cell = UnitCell(1,1,1)
@pytest.mark.parametrize("operations_string, symmetric, asymmetric", ops, sym, asym)
def test_exchange_matrix_constraints_object_gives_correct_results(operations_string, symmetric, asymmetric):
    operations = string_to_ops(operations_string)
    transforms = [op.point_operation_in_cartesian(cell) for op in operations]
    constraints = ExchangeMatrixConstraints(transforms)
    
    expressions = constraints._matrix_form_strings()

