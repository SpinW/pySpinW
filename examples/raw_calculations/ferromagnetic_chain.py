"""Basic Heisenberg ferromagnetic chain example.

See https://spinw.org/tutorials/01tutorial
"""
import sys

import numpy as np

from examples.raw_calculations.utils import run_example, plot, py_classes

def heisenberg_ferromagnet(n_q = 100, classes = py_classes):
    """Basic ferromagnet."""
    rust_kw = {'dtype':complex, 'order':'F'}
    Coupling = classes.coupling

    q_mags = np.linspace(0, 1, n_q).reshape(-1, 1)
    q_vectors = np.array([0, 1, 0]).reshape(1, 3) * q_mags
    positions = [np.array([0., 0., 0.])]  # single site at origin

    # Single site
    rotations = [np.eye(3, **rust_kw)]
    magnitudes = np.array([1.0])  # spin-1
    couplings = [Coupling(0, 0, np.eye(3, **rust_kw), inter_site_vector=np.array([0., 1., 0.])),
                 Coupling(0, 0, np.eye(3, **rust_kw), inter_site_vector=np.array([0., -1., 0.])),
                 ]

    return (rotations, magnitudes, q_vectors, couplings, positions)

if __name__ == "__main__":

    import matplotlib.pyplot as plt

    if len(sys.argv) > 1:
        use_rust = "py" not in sys.argv[1]
    else:
        use_rust = True

    q_mags = np.linspace(0, 1, 100)

    _, energies, sqw = run_example(heisenberg_ferromagnet, use_rust)

    # Note: we get complex data types with real part zero
    plot(q_mags, energies, sqw)

    # Compare with tutorial 1, 3rd last figure (https://spinw.org/tutorial1_05.png)
    plt.show()
