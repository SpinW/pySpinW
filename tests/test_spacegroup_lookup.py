import pytest
from pyspinw.interface import spacegroup
from pyspinw.symmetry.group import NoSuchGroup

expected_lookups = [
    ("c2/m", 12, "C 1 2/m 1", None),
    #("b2/m", 12, "B 1 1 2/m", None),
    ("P4/n:1", 85, "P 4/n", "1"),
    ("P4/n:2", 85, "P 4/n", "2"),
    ("P1", 1, "P 1", None)
]

@pytest.mark.parametrize("search_string, number, full_name, choice", expected_lookups)
def test_expected_lookups(search_string, number, full_name, choice):
    """ Check that the lookup works for some examples"""
    group = spacegroup(search_string)

    assert group.number == number
    assert group.symbol == full_name
    if choice is not None:
        assert group.choice == choice


# List of terms that should raise NoSuchGroup errors along with expected suggestions
expected_errors = [
    ("p7", []),
    ("", []),
    ("zzzzz", [])
]

@pytest.mark.parametrize("search_string, required_suggestions", expected_errors)
def test_expected_errors(search_string, required_suggestions):
    """ Check that ambiguous or wrong names throw errors, with appropriate suggestions in the error"""
    if expected_errors:
        for term in required_suggestions:
            with pytest.raises(NoSuchGroup, match=term):
                spacegroup(search_string)
    else:
        with pytest.raises(NoSuchGroup):
            spacegroup(search_string)
