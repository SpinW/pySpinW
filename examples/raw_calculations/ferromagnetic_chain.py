"""Basic Heisenberg ferromagnetic chain example.

See https://spinw.org/tutorials/01tutorial
"""
import sys

import numpy as np

try:
    from pyspinw.rust import spinwave_calculation as rs_spinwave, Coupling as RsCoupling
    RUST_AVAILABLE = True
except ModuleNotFoundError:
    RUST_AVAILABLE = False

from pyspinw.calculations.spinwave import spinwave_calculation as py_spinwave, Coupling as PyCoupling

def heisenberg_ferromagnet(n_q = 100, rust = False):
    """Basic ferromagnet."""
    if rust:
        rust_kw = {'dtype':complex, 'order':'F'}
        Coupling = RsCoupling
    else:
        rust_kw = {}
        Coupling = PyCoupling

    q_mags = np.linspace(0, 1, n_q).reshape(-1, 1)
    q_vectors = np.array([0, 1, 0]).reshape(1, 3) * q_mags

    # Single site
    rotations = [np.eye(3, **rust_kw)]
    magnitudes = np.array([1.0])  # spin-1
    couplings = [Coupling(0, 0, np.eye(3, **rust_kw), inter_site_vector=np.array([0., 1., 0.])),
                 Coupling(0, 0, np.eye(3, **rust_kw), inter_site_vector=np.array([0., -1., 0.])), ]

    return (rotations, magnitudes, q_vectors, couplings)

if __name__ == "__main__":

    import matplotlib.pyplot as plt

    n_q = 100
    use_rust = ("py" not in sys.argv[1]) if len(sys.argv) > 1 else RUST_AVAILABLE
    spinwave_calculation = rs_spinwave if use_rust else py_spinwave
    q_mags = np.linspace(0, 1, n_q)

    structure = heisenberg_ferromagnet(n_q, rust=use_rust)
    energies = spinwave_calculation(*structure)

    # Note: we get complex data types with real part zero

    plt.plot(q_mags, energies)

    # Compare with tutorial 1, 3rd last figure (https://spinw.org/tutorial1_05.png)
    plt.show()
