"""Example of a Kagome supercell.

This is the magnet in MATLAB SpinW tutorial 8:
    https://spinw.org/tutorials/08tutorial
"""
import sys

import numpy as np

from examples.raw_calculations.utils import run_example, plot, py_classes
from pyspinw.hamiltonian import omegasum


def kagome_supercell(n_q = 100, classes = py_classes):
    """A sqrt(3) x sqrt(3) Kagome antiferromagnet supercell lattice."""
    rust_kw = {'dtype':complex, 'order':'F'}
    Coupling = classes.coupling


    # we index the supercell by indexing each unit cell in order: so that the
    # 'central' atom is 0 mod 3, the top-left atom is 1 mod 3, the right is 2 mod 3
    # so we can create the unit cell couplings and then we just need to add
    # the remaining couplings between unit cells
    # here over the 3x3 supercell it is done as
    #
    # 6 \ 7 \ 8
    # ----------
    #  3 \ 4 \ 5
    #  ----------
    #   0 \ 1 \ 2
    #
    # Within each unit cell there are 3 spins, at (0.5, 0), (0, 0.5) and (0.5, 0.5)
    s3 = np.sqrt(3) / 2
    def down():
        # ↙️
        return np.array([[s3, 0, -0.5], [-0.5, 0, -s3], [0, 1, 0]], dtype=complex, order='F')

    def up():
        # ↖
        return np.array([[s3, 0, -0.5], [0.5, 0, s3], [0, -1, 0]], dtype=complex, order='F')

    def right():
        # →
        return np.array([[0, 0, 1.], [1., 0, 0], [0, 1, 0]], dtype=complex, order='F')

    rotations = [
            up(),    # Cell 0, (0.5, 0)
            up(),    # Cell 0, (0, 0.5)
            down(),  # Cell 0, (0.5, 0.5)
            right(), # Cell 1, (0.5, 0)
            right(), # Cell 1, (0, 0.5)
            up(),    # Cell 1, (0.5, 0.5)
            down(),
            down(),
            right(),
            right(),
            right(),
            up(),
            down(),
            down(),
            right(),
            up(),
            up(),
            down(),
            down(),
            down(),
            right(),
            up(),
            up(),
            down(),
            right(),
            right(),
            up(),
        ]
    magnitudes = np.array([1] * 27)  # spin 1
    unit_cell_positions = [np.array([1/6., 0., 0.]), np.array([0., 1/6., 0.]), np.array([1/6., 1/6., 0.])]
    positions = (unit_cell_positions  # unit cell 0
                + [pos + np.array([1., 0., 0.])/3. for pos in unit_cell_positions]  # unit cell 1
                + [pos + np.array([2., 0., 0.])/3. for pos in unit_cell_positions]  # unit cell 2
                + [pos + np.array([0., 1., 0.])/3. for pos in unit_cell_positions]  # unit cell 3
                + [pos + np.array([1., 1., 0.])/3. for pos in unit_cell_positions]  # unit cell 4
                + [pos + np.array([2., 1., 0.])/3. for pos in unit_cell_positions]  # unit cell 5
                + [pos + np.array([0., 2., 0.])/3. for pos in unit_cell_positions]  # unit cell 6
                + [pos + np.array([1., 2., 0.])/3. for pos in unit_cell_positions]  # unit cell 7
                + [pos + np.array([2., 2., 0.])/3. for pos in unit_cell_positions]) # unit cell 8

    # The inter_site_vector is actually the vector between unit cells (see eq 10 and eq 14 in Toth+Lake)
    # So we define the pairs within each cell first - there are 2 for each cell
    other_pairs = [
        (0, 2), (1, 2),      (3, 5), (4, 5),      (6, 8), (7, 8),     # Cells 0-2
        (9, 11), (10, 11),   (12, 14), (13, 14),  (15, 17), (16, 17), # Cells 3-5
        (18, 20), (19, 20),  (21, 23), (22, 23),  (24, 26), (25, 26), # Cells 6-8
        (1, 9), (2, 9),      (4, 12), (5, 12),    (7, 15), (8, 15),   # Cells 0-2 vertical supercell
        (10, 18), (11, 18),  (13, 21), (14, 21),  (16, 24), (17, 24), # Cells 3-5 vertical supercell
        (0, 4), (2, 4),      (3, 7), (5, 7),                          # Cells 0-2 horizontal supercell
        (9, 13), (11, 13),   (12, 16), (14, 16),                      # Cells 3-5 horizontal supercell
        (18, 22), (20, 22),  (21, 25), (23, 25),                      # Cells 6-8 horizontal supercell
    ]
    # Now, couplings between unit cells parallel to the a-axis (also two per cell)
    horizontal_pairs = [
        # Cells 0-2          Cells 3-5            Cells 6-8
        (6, 1), (8, 1),      (15, 10), (17, 10),  (24, 19), (26, 19),
    ]
    # Now, couplings between unit cells parallel to the b-axis (also two per cell)
    vertical_pairs = [
        (19, 0), (20, 0),    (22, 3), (23, 3),    (25, 6), (26, 6),   # Cells 6-8
    ]

    couplings = []

    couplings.extend(Coupling(idx1, idx2, np.eye(3, **rust_kw), np.array([1., 0, 0])) for (idx1, idx2) in horizontal_pairs)
    couplings.extend(Coupling(idx2, idx1, np.eye(3, **rust_kw), np.array([-1., 0, 0])) for (idx1, idx2) in horizontal_pairs)
    couplings.extend(Coupling(idx1, idx2, np.eye(3, **rust_kw), np.array([0., 1, 0])) for (idx1, idx2) in vertical_pairs)
    couplings.extend(Coupling(idx2, idx1, np.eye(3, **rust_kw), np.array([0., -1, 0])) for (idx1, idx2) in vertical_pairs)
    couplings.extend(Coupling(idx1, idx2, np.eye(3, **rust_kw), np.array([0., 0, 0])) for (idx1, idx2) in other_pairs)
    couplings.extend(Coupling(idx2, idx1, np.eye(3, **rust_kw), np.array([0., 0, 0])) for (idx1, idx2) in other_pairs)

    q_mags = 0.5 * np.linspace(0, 1, n_q + 1).reshape(-1, 1)

    q_vectors = np.concatenate(
        (
            q_mags[::-1].reshape(-1, 1) * np.array([-3, 0, 0]).reshape(1, -1),
            q_mags[1:].reshape(-1, 1) * np.array([3, 3, 0]).reshape(1, -1),
        )
    )

    rlu2cart = np.sqrt([[1./36, 1./108, 0], [0, 1./27, 0], [0, 0, 1./1600]]) * 2 * np.pi
    return rotations, magnitudes, q_vectors, couplings, positions, rlu2cart

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    if len(sys.argv) > 1:
        use_rust = "py" not in sys.argv[1]
    else:
        use_rust = True

    structure, energies, sqw = run_example(kagome_supercell, use_rust)

    q_vectors = structure[2]
    indices = np.arange(201)
    label_indices = [0, 100, 200]
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
    fg.axes[1].set_ylim(0, 1)
    plt.show()
