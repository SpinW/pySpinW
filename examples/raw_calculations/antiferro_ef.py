"""Example of an antiferromagnetic chain with an external magnetic field.

Based on the (currently unpublished) spinW Tutorial 29:
https://github.com/SpinW/spinw/blob/dc06a4bec2c44bcfde6baa630ecb6329fb4b68ba/tutorials/tutorial/tutorial29.m
"""
import sys

import numpy as np
from pyspinw.calculations.spinwave import Coupling as PyCoupling, MagneticField as PyField

from examples.raw_calculations.utils import run_example

def antiferro_ef(n_q = 100, coupling_class = PyCoupling, field_class = PyField):
    """Antiferromagnetic chain.

    We use a 2x1x1 supercell to capture the magnetic rotation periodicity.
    """
    rust_kw = {'dtype':complex, 'order':'F'}
    Coupling = coupling_class
    MagneticField = field_class

    rotations = [np.eye(3, **rust_kw), np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]], **rust_kw)]
    magnitudes = np.array([1.0]*2)

    aniso_array = np.diag(np.array([0., 0., -0.1], **rust_kw))

    couplings = [
        # bonds
        Coupling(0, 1, np.eye(3, **rust_kw), inter_site_vector=np.array([0., 1., 0.])),
        Coupling(0, 1, np.eye(3, **rust_kw), inter_site_vector=np.array([0., -1., 0.])),
        Coupling(1, 0, np.eye(3, **rust_kw), inter_site_vector=np.array([0., 1., 0.])),
        Coupling(1, 0, np.eye(3, **rust_kw), inter_site_vector=np.array([0., -1., 0.])),
        # anisotropies
        Coupling(0, 0, aniso_array, inter_site_vector=np.array([0., 0., 0.])),
        Coupling(1, 1, aniso_array, inter_site_vector=np.array([0., 0., 0.])),
    ]

    g_tensors = [np.eye(3, **rust_kw) * 2] * 2

    ext_field = MagneticField(vector=np.array([0., 0., 7.], **rust_kw),
                              g_tensors=g_tensors)

    q_mags = np.linspace(0, 1, n_q).reshape(-1, 1)
    q_vectors = np.array([0, 1, 0]).reshape(1, 3) * q_mags

    return (rotations, magnitudes, q_vectors, couplings, ext_field)

if __name__ == "__main__":

    import matplotlib.pyplot as plt

    if len(sys.argv) > 1:
        use_rust = "py" not in sys.argv[1]
    else:
        use_rust = True

    _, energies = run_example(antiferro_ef, use_rust, has_field=True)


    # Note: we get complex data types with real part zero

    plt.plot(np.linspace(0, 1, 100).reshape(-1, 1), energies)

    #plt.savefig("fig.png")
    # Compare with tutorial 2, second last figure (https://spinw.org/tutorial2_05.png)
    plt.show()
