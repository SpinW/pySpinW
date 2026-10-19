import pytest

from pyspinw import LatticeSite
from pyspinw.anisotropy import Anisotropy, _all_anisotropies

import numpy as np

matrices = [
    np.array([[1,2,3],
             [2,1,2],
             [3,2,1]]), # General

    np.diag([1,1,-1]), # General

    np.diag([1,0,0]) # Angle-magnitude
]

sites = [
    LatticeSite(1/2, 1/2, 1/2, 0, 1, 0),
    LatticeSite(0,0,0,1,0,1)

]

def test_all_classes_present():
    """ Test that checks that we are testing all the anisotropy classes when generating using matrix -> specialise"""

    instances = [Anisotropy(sites[0], anisotropy_matrix=matrix).specialise() for matrix in matrices]

    for instance in instances:
        print(instance)

    classes = [instance.__class__ for instance in instances]

    for target_class in _all_anisotropies:
        assert target_class in classes

@pytest.mark.parametrize("name", ["Bob", "Dave"])
@pytest.mark.parametrize("site", sites)
@pytest.mark.parametrize("matrix", matrices)
def test_anisotropy_serialisation(site, matrix, name):
    anisotropy = Anisotropy(site, anisotropy_matrix=matrix, name=name).specialise()

    json = anisotropy.serialise()
    deserialised = Anisotropy.deserialise(json)

    assert anisotropy.__class__ == deserialised.__class__

    assert np.all(anisotropy._anisotropy_matrix == deserialised._anisotropy_matrix)
    assert anisotropy.name == deserialised.name