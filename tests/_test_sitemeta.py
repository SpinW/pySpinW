import pytest

from pyspinw.sitemeta import SiteMetadata
from ase.data import chemical_symbols

@pytest.mark.parametrize("rest", "agAGxX_-[] ")
@pytest.mark.parametrize("element", chemical_symbols[1:])
def test_element_assignment_element_exists(element, rest):
    test_string = element + rest
    if test_string not in chemical_symbols:
        meta = SiteMetadata.metadata_from_name_and_spin(test_string)
        assert meta.element == element


@pytest.mark.parametrize("chr2", "XYZ")
@pytest.mark.parametrize("chr1", ["", " ", "x", "y", "z"])
def test_element_assignment_no_element(chr1, chr2):
    test_string = chr1 + chr2
    if test_string.strip() not in chemical_symbols:
        assert SiteMetadata.metadata_from_name_and_spin(test_string).element is None
