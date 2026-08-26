""" Test for the operation that converts desciptions of spacegroup operations (e.g. "z,-y,x") into objects """

import pytest

from pyspinw.symmetry.group import database
from pyspinw.symmetry.operations import SpaceOperation

operations = set([op for sg in database.spacegroups for op in sg.operations])

@pytest.mark.parametrize("op", operations, ids=lambda x: f"[{x.text_form}]")
def test_operation_from_text(op):
    """ Check that parsing the text form returns the same operation (or at least, one numerically the same)"""

    new_op = SpaceOperation.from_text(op.text_form)

    assert new_op == op, "Operators should be identical"