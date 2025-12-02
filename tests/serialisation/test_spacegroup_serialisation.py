""" Tests for the serialisation of spacegroups """

import pytest

from pyspinw.symmetry.group import database, SpaceGroup

@pytest.mark.parametrize("group", database.spacegroups)
def test_serialisation_string(group):
    """ Check that the serialisation string is something that we can use to lookup spacegroups correctly with """
    string = group._serialisation_string()
    from_string = database.spacegroup_by_name(string)

    assert from_string.symbol == group.symbol
    assert from_string.choice == group.choice
    assert from_string.short_symbol == group.short_symbol


