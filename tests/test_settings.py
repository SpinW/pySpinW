import itertools

import pytest

from pyspinw.symmetry.group import database
from pyspinw.symmetry.settings import UniqueAxis, AxisPermutation, Setting


@pytest.mark.parametrize("string", ["a", "b", "c", "-a", "-b", "-c"])
def test_unique_axis(string):
    test_string = UniqueAxis.from_string(string).to_string()

    assert test_string == string


# create all permuations of "abc" and "ab-c"
perms = list(itertools.permutations("abc"))
perms += [["-c" if s == "c" else s for s in perm] for perm in perms]
perms = ["".join(perm) for perm in perms]

@pytest.mark.parametrize("string", perms)
def test_permutation(string):
    """ Test that perutations convert to and from strings consistently"""
    test_string = AxisPermutation.from_string(string).to_string()

    assert test_string == string

@pytest.mark.parametrize("group", database.spacegroups)
def test_choices(group):
    """ Test that the whole settings object goes from and to strings consistently """
    string = "" if group.choice is None else group.choice
    test_string = Setting.from_string(string).to_string()
    assert string == test_string
