import numpy as np

from pyspinw.coupling import HeisenbergCoupling, DMCoupling
from pyspinw.site import LatticeSite

a = LatticeSite(0,0,0, name="A")
b = LatticeSite(0.5,0.5,0.5, name="B")

def test_create_heisenberg():
    hc = HeisenbergCoupling(site_1=a, site_2=b, j=10, name="J1")

    assert np.all(np.abs(hc.coupling_matrix - 10*np.eye(3)) < 1e-10)


def test_create_dm():
    hc = DMCoupling(site_1=a, site_2=b, d_x=10, d_y=4, d_z=1, name="D1")

    expected = np.array([
        [0, 1, -4],
        [-1, 0, 10],
        [4, -10, 0]
    ])

    assert np.all(np.abs(hc.coupling_matrix - expected) < 1e-10)

