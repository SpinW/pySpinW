import pytest

from pyspinw.symmetry.group import database
from pyspinw.symmetry.operations import SpaceOperation

operations = set([op for sg in database.spacegroups for op in sg.operations])

@pytest.mark.parametrize("op", operations, ids=lambda x: f"[{x.text_form}]")
def test_operation_from_text(op):
    new_op = SpaceOperation.from_text(op.text_form)

    assert new_op == op, "Operators should be identical"