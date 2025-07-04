"""Example of an antiferromagnetic chain.

See https://spinw.org/tutorials/02tutorial
"""
import numpy as np

try:
    from pyspinw.rust import spinwave_calculation, Coupling
except ModuleNotFoundError:
    from pyspinw.calculations.spinwave import spinwave_calculation, Coupling

def antiferro_chain(n_q = 100):
    """Antiferromagnetic chain.

    We use a 2x1x1 supercell to capture the magnetic rotation periodicity.
    """
    rust_kw = {'dtype':complex, 'order':'F'}
    rotations = [np.eye(3, **rust_kw), np.array([[-1, 0, 0], [0, 1, 0], [0, 0, -1]], **rust_kw)]
    magnitudes = np.array([1.0]*2)

    couplings = [
        Coupling(0, 1, np.eye(3, **rust_kw), inter_site_vector=np.array([0., 1., 0.])),
        Coupling(0, 1, np.eye(3, **rust_kw), inter_site_vector=np.array([0., -1., 0.])),
        Coupling(1, 0, np.eye(3, **rust_kw), inter_site_vector=np.array([0., 1., 0.])),
        Coupling(1, 0, np.eye(3, **rust_kw), inter_site_vector=np.array([0., -1., 0.])),
    ]

    q_mags = np.linspace(0, 1, n_q).reshape(-1, 1)
    q_vectors = np.array([0, 1, 0]).reshape(1, 3) * q_mags

    return (rotations, magnitudes, q_vectors, couplings)

if __name__ == "__main__":

    import matplotlib.pyplot as plt

    structure = antiferro_chain()
    energies = spinwave_calculation(*structure)

    # Note: we get complex data types with real part zero

    plt.plot(np.linspace(0, 1, 100).reshape(-1, 1), energies)

    #plt.savefig("fig.png")
    # Compare with tutorial 2, second last figure (https://spinw.org/tutorial2_05.png)
    plt.show()
