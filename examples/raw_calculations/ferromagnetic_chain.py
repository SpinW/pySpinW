"""Basic Heisenberg ferromagnetic chain example.

See https://spinw.org/tutorials/01tutorial
"""
import numpy as np

try:
    from pyspinw.rust import Coupling, spinwave_calculation
except ModuleNotFoundError:
    from pyspinw.calculations.spinwave import Coupling, spinwave_calculation

def heisenberg_ferromagnet(n_q = 100):
    """Basic ferromagnet."""
    q_mags = np.linspace(0, 1, n_q).reshape(-1, 1)
    q_vectors = np.array([0, 1, 0]).reshape(1, 3) * q_mags

    # Single site
    rotations = [np.eye(3, dtype=complex, order='F')]
    magnitudes = np.array([1.0])  # spin-1
    couplings = [Coupling(0, 0, np.eye(3, dtype=complex, order='F'), inter_site_vector=np.array([0., 1., 0.])),
                 Coupling(0, 0, np.eye(3, dtype=complex, order='F'), inter_site_vector=np.array([0., -1., 0.])), ]

    return (rotations, magnitudes, q_vectors, couplings)

if __name__ == "__main__":

    import matplotlib.pyplot as plt

    n_q = 100
    structure = heisenberg_ferromagnet(n_q)
    q_mags = np.linspace(0, 1, n_q)
    energies = spinwave_calculation(*structure)

    # Note: we get complex data types with real part zero

    plt.plot(q_mags, energies)

    # Compare with tutorial 1, 3rd last figure (https://spinw.org/tutorial1_05.png)
    plt.show()
