import pytest
from pyspinw import SiteMetadata, LatticeSite
import numpy as np

element_list = [None, "C", "B", "Na", "Fe"]
name_list = [None, "A", "B", "C", "Na", "X", "Y", "Z", "Ca 2+", "Cl-"]
color_list = [None, (1,1,1), (0.5, 0.1, 0.0)]
radius_list = [None, 1, 0.5, 10.7895641361]

@pytest.mark.parametrize("radius", radius_list)
@pytest.mark.parametrize("color", color_list)
@pytest.mark.parametrize("element", element_list)
def test_site_metadata_serialisation(element, color, radius):
    """ Check that metadata objects serialise correctly """
    meta = SiteMetadata(element=element, radius=radius, color=color)
    json = meta.serialise()
    deserialised = SiteMetadata.deserialise(json)

    assert isinstance(deserialised, SiteMetadata)

    assert meta.element == deserialised.element
    assert np.all(np.array(meta.color) == np.array(deserialised.color))
    assert meta.radius == deserialised.radius


@pytest.mark.parametrize("name", name_list)
def test_site_serialisation_with_implied_metadata(name):
    """ Check site serialisation with metadata implied by name """
    site = LatticeSite(1/2, 1/2, 1/2, name=name)
    json = site.serialise()
    deserialised = LatticeSite.deserialise(json)

    assert isinstance(deserialised, LatticeSite)

    assert site.name == deserialised.name
    assert site.metadata.element == deserialised.metadata.element
    assert np.all(np.array(site.metadata.color) == np.array(deserialised.metadata.color))
    assert site.metadata.radius == deserialised.metadata.radius

@pytest.mark.parametrize("name", name_list)
@pytest.mark.parametrize("color", color_list)
def test_site_serialisation_with_color(name, color):
    """ Check site serialisation with metadata implied by name and color keyword """
    site = LatticeSite(1/2, 1/2, 1/2, name=name, color=color)
    json = site.serialise()
    deserialised = LatticeSite.deserialise(json)

    assert isinstance(deserialised, LatticeSite)

    assert site.name == deserialised.name
    assert site.metadata.element == deserialised.metadata.element
    assert np.all(np.array(site.metadata.color) == np.array(deserialised.metadata.color))
    assert site.metadata.radius == deserialised.metadata.radius


@pytest.mark.parametrize("radius", radius_list)
@pytest.mark.parametrize("color", color_list)
@pytest.mark.parametrize("element", element_list)
@pytest.mark.parametrize("name", name_list)
def test_site_serialisation_with_metadata(name, element, color, radius):
    """ Check site serialisation with metadata implied by name and metadata keyword """

    meta = SiteMetadata(element=element, radius=radius, color=color)

    site = LatticeSite(1/2, 1/2, 1/2, name=name, metadata=meta)
    json = site.serialise()
    deserialised = LatticeSite.deserialise(json)

    assert isinstance(deserialised, LatticeSite)

    assert site.name == deserialised.name
    assert site.metadata.element == deserialised.metadata.element
    assert np.all(np.array(site.metadata.color) == np.array(deserialised.metadata.color))
    assert site.metadata.radius == deserialised.metadata.radius