""" Check that the anisotropy specialisation methods work """

import pytest
import numpy as np

from pyspinw import AxisMagnitudeAnisotropy, LatticeSite
from pyspinw.anisotropy import specialise_anisotropy

n=10
filters = [
    [1, 1, 1],
    [1, 1 ,0],
    [1, 0 ,1],
    [0, 1 ,1],
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1],
]

rng = np.random.default_rng(668)

@pytest.mark.parametrize("filter", filters, ids=lambda x: f"filter={repr(x)}")
def test_axis_magnitude_specialisation_recovers(filter):
    """ Check that axis-magnitude anisotropies """

    filter = np.array(filter, dtype=float)
    site = LatticeSite(0,0,0)

    for i in range(n):
        vec = rng.normal(0,1, (3,))
        vec /= np.sqrt(np.sum(vec**2))

        vec *= filter

        a = 0
        while a == 0:
            a = rng.normal(0,1)


        original = AxisMagnitudeAnisotropy(site=site, a=a, direction=vec)

        rebuilt = specialise_anisotropy(original.generalise())

        assert isinstance(rebuilt, AxisMagnitudeAnisotropy), ("Anisotropy should be converted back "
                                                              "to AxisMagnitudeAnisotropy.")

        assert (np.allclose(original.direction, rebuilt.direction) or
                np.allclose(original.direction, -rebuilt.direction)), ("Recovered direction should be the same as"
                                                                       "the original, up to sign.")

        assert np.allclose(original.a, rebuilt.a), "Recovered `a` constant should be the same as the original."



