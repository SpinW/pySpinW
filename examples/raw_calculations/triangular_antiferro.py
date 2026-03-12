"""Example of a triangular antiferromagnet using the rotating frame method

This is the magnet in MATLAB SpinW tutorial 12:
    https://spinw.org/tutorials/12tutorial
"""
import sys

import numpy as np

from examples.raw_calculations.utils import run_example, plot, py_classes
from pyspinw.hamiltonian import omegasum


def triangular_antiferro(n_q = 100, classes = py_classes):
    """A sqrt(3) x sqrt(3) Kagome antiferromagnet supercell lattice."""
    rust_kw = {'dtype':complex, 'order':'F'}
    Coupling = classes.coupling

    rotations = [np.array([[1, 0, 0], [0, 0, 1], [0, 1, 0]], **rust_kw)]
    magnitudes = np.array([1.5])
    positions = [np.array([0., 0., 0.]),]
    couplings = [Coupling(0, 0, np.eye(3, **rust_kw), inter_site_vector=np.array([ 0.,  1., 0.])),
                 Coupling(0, 0, np.eye(3, **rust_kw), inter_site_vector=np.array([ 0., -1., 0.])),
                 Coupling(0, 0, np.eye(3, **rust_kw), inter_site_vector=np.array([ 1.,  0., 0.])),
                 Coupling(0, 0, np.eye(3, **rust_kw), inter_site_vector=np.array([-1.,  0., 0.])),
                 Coupling(0, 0, np.eye(3, **rust_kw), inter_site_vector=np.array([ 1.,  1., 0.])),
                 Coupling(0, 0, np.eye(3, **rust_kw), inter_site_vector=np.array([-1., -1., 0.])),
                 # Final "coupling" is a planar single-ion anisotropy
                 Coupling(0, 0, np.diag(np.array([0, 0, 0.1], **rust_kw)), inter_site_vector=np.array([0., 0., 0.]))]
    rlu_to_cart = np.array([[2.0943951, 1.20919958, 0.], [0., 2.41839915, 0.], [0., 0., 1.57079633]], order='F')
    field = None
    rotating_frame = [np.array([1./3, 1./3, 1.]), np.array([0., 0., 1.])]

    q_mags = np.linspace(0.01, 0.99, n_q + 1).reshape(-1, 1)
    q_vectors = q_mags.reshape(-1, 1) * np.array([1, 1, 0]).reshape(1, -1)

    return rotations, magnitudes, q_vectors, couplings, positions, rlu_to_cart, field, rotating_frame

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    if len(sys.argv) > 1:
        use_rust = "py" not in sys.argv[1]
    else:
        use_rust = True

    structure, energies, sqw = run_example(triangular_antiferro, use_rust)

    q_vectors = structure[2]
    indices = np.arange(101)
    label_indices = [0, 100]
    labels = [str(q_vectors[idx, :]) for idx in label_indices]

    energies = [np.sort(energy.real) for energy in energies]
    positive_energies = [energy[energy > 0] for energy in energies]
    min_energy = min([np.min(energy) for energy in positive_energies])
    translated_energies = [energy - min_energy for energy in positive_energies]

    # Note: we get complex data types with real part zero
    translated_energies, sqw = omegasum(translated_energies, sqw)

    fg = plot(indices, translated_energies, sqw, show=False)
    # plt.plot(indices, [method.value for method in result.method])

    #plt.savefig("fig.png")
    fg.axes[1].set_ylim(0, 6)
    plt.show()
