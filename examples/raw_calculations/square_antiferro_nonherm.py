"""Example of a square lattice antiferromagnet with axial SIA and imaginary modes
   
   This example is to test the non-hermitian calculations
"""
import sys

import numpy as np

from examples.raw_calculations.utils import run_example, py_classes


def square_antiferro_nonherm(n_q = 100, classes = py_classes):
    """A square lattice antiferromagnet in the a-b plane with large axial SIA || c"""
    rust_kw = {'dtype':complex, 'order':'F'}
    Coupling = classes.coupling

    rotations = [np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], **rust_kw),
                 np.array([[0, 0, -1], [1, 0, 0], [0, -1, 0]], **rust_kw),
                 np.array([[0, 0, -1], [1, 0, 0], [0, -1, 0]], **rust_kw),
                 np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], **rust_kw)]
    magnitudes = np.array([0.5] * 4)
    positions = [np.array([0., 0., 0.]), np.array([0.5, 0, 0]), np.array([0, 0.5, 0]), np.array([0.5, 0.5, 0])]
    couplings = [Coupling(0, 1, np.eye(3, **rust_kw), inter_site_vector=np.array([ 0.,  0., 0.])),
                 Coupling(0, 2, np.eye(3, **rust_kw), inter_site_vector=np.array([ 0.,  0., 0.])),
                 Coupling(1, 0, np.eye(3, **rust_kw), inter_site_vector=np.array([ 1.,  0., 0.])),
                 Coupling(1, 3, np.eye(3, **rust_kw), inter_site_vector=np.array([ 0.,  0., 0.])),
                 Coupling(2, 3, np.eye(3, **rust_kw), inter_site_vector=np.array([ 0.,  0., 0.])),
                 Coupling(2, 0, np.eye(3, **rust_kw), inter_site_vector=np.array([ 0.,  1., 0.])),
                 Coupling(3, 2, np.eye(3, **rust_kw), inter_site_vector=np.array([ 1.,  0., 0.])),
                 Coupling(3, 1, np.eye(3, **rust_kw), inter_site_vector=np.array([ 0.,  1., 0.])),
                 Coupling(1, 0, np.eye(3, **rust_kw), inter_site_vector=np.array([ 0.,  0., 0.])),
                 Coupling(2, 0, np.eye(3, **rust_kw), inter_site_vector=np.array([ 0.,  0., 0.])),
                 Coupling(0, 1, np.eye(3, **rust_kw), inter_site_vector=np.array([-1.,  0., 0.])),
                 Coupling(3, 1, np.eye(3, **rust_kw), inter_site_vector=np.array([ 0.,  0., 0.])),
                 Coupling(3, 2, np.eye(3, **rust_kw), inter_site_vector=np.array([ 0.,  0., 0.])),
                 Coupling(0, 2, np.eye(3, **rust_kw), inter_site_vector=np.array([ 0., -1., 0.])),
                 Coupling(2, 3, np.eye(3, **rust_kw), inter_site_vector=np.array([-1.,  0., 0.])),
                 Coupling(1, 3, np.eye(3, **rust_kw), inter_site_vector=np.array([ 0., -1., 0.])),
                 # Final set of "couplings" are the axial single-ion anisotropy
                 Coupling(0, 0, np.diag(np.array([0, 0, -0.4], **rust_kw)), inter_site_vector=np.array([0., 0., 0.])),
                 Coupling(1, 1, np.diag(np.array([0, 0, -0.4], **rust_kw)), inter_site_vector=np.array([0., 0., 0.])),
                 Coupling(2, 2, np.diag(np.array([0, 0, -0.4], **rust_kw)), inter_site_vector=np.array([0., 0., 0.])),
                 Coupling(3, 3, np.diag(np.array([0, 0, -0.4], **rust_kw)), inter_site_vector=np.array([0., 0., 0.]))]

    q_mags = np.linspace(0., 2., n_q + 1).reshape(-1, 1)
    q_vectors = q_mags.reshape(-1, 1) * np.array([1, 1, 0]).reshape(1, -1)

    return rotations, magnitudes, q_vectors, couplings, positions

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    if len(sys.argv) > 1:
        use_rust = "py" not in sys.argv[1]
    else:
        use_rust = True

    structure, energies, sqw = run_example(square_antiferro_nonherm, use_rust)

    q_vectors = structure[2]
    indices = np.arange(101)
    label_indices = [0, 100]
    labels = [str(q_vectors[idx, :]) for idx in label_indices]

    fg = plt.figure(figsize=(8, 6))
    plt.subplot(2, 1, 1)
    plt.plot(indices, np.real(energies), '.k')
    plt.plot(indices, np.imag(energies), 'or')
    plt.xlabel("Wavevector (r.l.u.)")
    plt.ylabel("Energy (meV)")

    plt.subplot(2, 1, 2)
    plt.plot(indices, sqw)
    plt.xlabel("Wavevector (r.l.u.)")
    plt.ylabel("S(q,ω) (arb. units)")

    plt.tight_layout()
    fg.axes[1].set_ylim(0, 10)
    plt.show()
