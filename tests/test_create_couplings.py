import numpy as np

from pyspinw.coupling import HeisenbergCoupling, DMCoupling


def test_create_heisenberg():
    hc = HeisenbergCoupling(site_1="A", site_2="B", j=10)

    assert np.all(np.abs(hc.coupling_matrix - 10*np.eye(3)) < 1e-10)


def test_create_dm():
    hc = DMCoupling(site_1="A", site_2="B", d_x=10, d_y=4, d_z=1)

    expected = np.array([
        [0, 1, -4],
        [-1, 0, 10],
        [4, -10, 0]
    ])

    assert np.all(np.abs(hc.coupling_matrix - expected) < 1e-10)

