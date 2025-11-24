""" Serialisation tests for sites """

import pytest

from pyspinw.site import LatticeSite, ImpliedLatticeSite

def test_serialisation():
    """ Check that the serialisation for lattice sites works with a few examples"""
    s1 = LatticeSite(1,2,3, 4, 5, 6, name="s1")
    s2 = LatticeSite(3,4,5, supercell_moments=[[7,8,9],[10,11,12]], name="s2")
    s3 = ImpliedLatticeSite(s1, 4,5,6,7,8,9)

    # Even though this doesn't make sense, it should serialise correctly
    s4 = ImpliedLatticeSite(s2, 2,3,4, supercell_moments=[1,2,3])

    for site in [s1, s2, s3, s4]:
        deserialised = LatticeSite.deserialise(site.serialise())

        assert site.i == deserialised.i
        assert site.j == deserialised.j
        assert site.k == deserialised.k

        assert site._moment_data == pytest.approx(deserialised._moment_data)
        assert site.name == deserialised.name

    for site, parent in [(s3, s1), (s4, s2)]:
        deserialised = LatticeSite.deserialise(site.serialise())

        assert site.parent_site.name == deserialised.parent_site.name
        assert parent.name == deserialised.parent_site.name

        assert parent.i == deserialised.parent_site.i
        assert parent.j == deserialised.parent_site.j
        assert parent.k == deserialised.parent_site.k

        assert parent._moment_data == pytest.approx(deserialised.parent_site._moment_data)
        assert parent.name == deserialised.parent_site.name
