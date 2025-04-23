import numpy as np
from pyspinw.calculations.spinwave import spinwave_calculation, Coupling

def heisenberg_ferromagnet():

    """
    Basic ferromagnet
    """

    q_mags = np.linspace(0, 1, 100).reshape(-1, 1)
    q_vectors = np.array([0, 1, 0]).reshape(1, 3) * q_mags

    # Single site
    rotations = np.eye(3).reshape(1, 3, 3)
    magnitudes = np.array([1.5])  # spin 3/2
    couplings = [Coupling(0, 0, np.eye(3), inter_site_vector=np.array([0, 1, 0])),
                 Coupling(0, 0, np.eye(3), inter_site_vector=np.array([0, -1, 0])),
                 ]

    energies = spinwave_calculation(rotations,
                                    magnitudes,
                                    q_vectors,
                                    couplings)

    return q_mags, energies

if __name__ == "__main__":

    import matplotlib.pyplot as plt

    q_mags, result = heisenberg_ferromagnet()

    # Note: we get complex data types with real part zero

    plt.plot(q_mags, result.raw_energies)

    plt.show()