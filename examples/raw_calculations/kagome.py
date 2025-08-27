"""Kagome ferromagnet example.

See https://spinw.org/tutorials/05tutorial.
"""
import sys

import numpy as np

try:
    from pyspinw.rust import spinwave_calculation as rs_spinwave, Coupling as RsCoupling
    RUST_AVAILABLE = True
except ModuleNotFoundError:
    RUST_AVAILABLE = False

from pyspinw.calculations.spinwave import spinwave_calculation as py_spinwave, Coupling as PyCoupling

def kagome_ferromagnet(n_q = 100, rust = False):
    """Basic ferromagnet on a kagome lattice."""
    if rust:
        rust_kw = {'dtype':complex, 'order':'F'}
        Coupling = RsCoupling
    else:
        rust_kw = {}
        Coupling = PyCoupling

    # Three sites, otherwise identical
    rotations = [np.eye(3, **rust_kw) for _ in range(3)]
    magnitudes = np.array([1.0]*3)  # spin-1

    # Each site coupled to two of each of the others, so there are 6 couplings,
    # And we need to add the other direction, so 12 entries in total
    #
    # 0-1 couplings have direction (0.5, 0,   0)
    # 0-2 couplings have direction (0.5, 0.5, 0)
    # 1-2 couplings have direction (0,   0.5, 0)
    #
    # We need to have them in both directions, for both permutations
    #  making 12 in total

    k = 0.5

    couplings = [
                 Coupling(0, 1, np.eye(3, **rust_kw), inter_site_vector=k*np.array([ 1.,  0., 0.])),
                 Coupling(0, 1, np.eye(3, **rust_kw), inter_site_vector=k*np.array([-1.,  0., 0.])),
                 Coupling(1, 0, np.eye(3, **rust_kw), inter_site_vector=k*np.array([ 1.,  0., 0.])),
                 Coupling(1, 0, np.eye(3, **rust_kw), inter_site_vector=k*np.array([-1.,  0., 0.])),
                 Coupling(0, 2, np.eye(3, **rust_kw), inter_site_vector=k*np.array([ 1.,  1., 0.])),
                 Coupling(0, 2, np.eye(3, **rust_kw), inter_site_vector=k*np.array([-1., -1., 0.])),
                 Coupling(2, 0, np.eye(3, **rust_kw), inter_site_vector=k*np.array([ 1.,  1., 0.])),
                 Coupling(2, 0, np.eye(3, **rust_kw), inter_site_vector=k*np.array([-1., -1., 0.])),
                 Coupling(1, 2, np.eye(3, **rust_kw), inter_site_vector=k*np.array([ 0.,  1., 0.])),
                 Coupling(1, 2, np.eye(3, **rust_kw), inter_site_vector=k*np.array([ 0., -1., 0.])),
                 Coupling(2, 1, np.eye(3, **rust_kw), inter_site_vector=k*np.array([ 0.,  1., 0.])),
                 Coupling(2, 1, np.eye(3, **rust_kw), inter_site_vector=k*np.array([ 0., -1., 0.]))
                 ]


    q_mags = 0.5*np.linspace(0, 1, n_q + 1).reshape(-1, 1)

    # q_vectors = np.concatenate((
    #         q_mags[::-1].reshape(-1, 1) * np.array([1, 0, 1]).reshape(1, -1),
    #         q_mags[1:].reshape(-1, 1) * np.array([0, 0, 1]).reshape(1, -1)
    # ))
    q_vectors = np.concatenate((
            q_mags[::-1].reshape(-1, 1) * np.array([-1, 0, 0]).reshape(1, -1),
            q_mags[1:].reshape(-1, 1) * np.array([1, 1, 0]).reshape(1, -1)
    ))

    # q_vectors = q_mags.reshape(-1, 1) * np.array([1, 1, 0]).reshape(1, -1)
    return rotations, magnitudes, q_vectors, couplings


if __name__ == "__main__":

    import matplotlib.pyplot as plt

    indices = np.arange(201)

    label_indices = [0, 100, 200]

    use_rust = ("py" not in sys.argv[1]) if len(sys.argv) > 1 else RUST_AVAILABLE
    spinwave_calculation = rs_spinwave if use_rust else py_spinwave

    rotations, magnitudes, q_vectors, couplings = kagome_ferromagnet(100, rust=use_rust)
    labels = [str(q_vectors[idx,:]) for idx in label_indices]
    energies = spinwave_calculation(rotations, magnitudes, q_vectors, couplings)

    plt.plot(indices, np.real(energies))
    # plt.plot(indices, [method.value for method in result.method])
    plt.xticks(label_indices, labels)

    # Compare with tutorial 5, 3rd last figure (https://spinw.org/tutorial5_05.png)
    plt.show()
