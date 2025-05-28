import numpy as np
from pyspinw.calculations.spinwave import spinwave_calculation, Coupling

def heisenberg_ferromagnet(n_q = 100):
    """
    Basic ferromagnet
    """

    q_mags = np.linspace(0, 1, n_q).reshape(-1, 1)
    q_vectors = np.array([0, 1, 0]).reshape(1, 3) * q_mags

    # Single site
    rotations = np.eye(3).reshape(1, 3, 3)
    magnitudes = np.array([1.5])  # spin 3/2
    couplings = [Coupling(0, 0, np.eye(3), inter_site_vector=np.array([0, 1, 0])),
                 Coupling(0, 0, np.eye(3), inter_site_vector=np.array([0, -1, 0])),
                 ]

    return (rotations, magnitudes, q_vectors, couplings)

if __name__ == "__main__":

    import matplotlib.pyplot as plt

    structure = heisenberg_ferromagnet()
    q_mags = np.linspace(0, 1, len(q_vectors))
    energies = spinwave_calculation(*structure)

    # Note: we get complex data types with real part zero

    plt.plot(q_mags, result.raw_energies)

    plt.show()
