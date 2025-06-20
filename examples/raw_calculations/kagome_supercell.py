"""Example of a Kagome supercell.

This is the magnet in MATLAB SpinW tutorial 8:
    https://spinw.org/tutorials/08tutorial
"""
import numpy as np

try:
    from pyspinw.rust import spinwave_calculation, Coupling
except ModuleNotFoundError:
    from pyspinw.calculations.spinwave import spinwave_calculation, Coupling


# define our rotation matrices
def rotation(theta):
    """Return the rotation matrix for an angle `theta` in the x-z plane."""
    return np.array(
        [
            [np.cos(theta), 0, np.sin(theta)],
            [np.sin(theta), 0, -np.cos(theta)],
            [0, 1, 0],
        ],
        dtype=complex, order='F',
    )


def kagome_supercell(n_q = 100):
    """A sqrt(3) x sqrt(3) Kagome antiferromagnet supercell lattice."""

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
    def down():
        # ↙️
        return rotation(-np.pi / 6)

    def up():
        # ↖
        return rotation(-5 * np.pi / 6)

    def right():
        # →
        return rotation(np.pi / 2)

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

    # The inter_site_vector is actually the vector between unit cells (see eq 10 and eq 14 in Toth+Lake)
    # So we define the pairs within each cell first - there are 2 for each cell
    other_pairs = [
        (0, 2), (1, 2),      (3, 5), (4, 5),      (6, 8), (7, 8),     # Cells 0-2
        (9, 11), (10, 11),   (12, 14), (13, 14),  (15, 17), (16, 17), # Cells 3-5
        (18, 20), (19, 20),  (21, 23), (22, 23),  (24, 26), (25, 26), # Cells 6-8
    ]
    # Now, couplings between unit cells parallel to the a-axis (also two per cell)
    horizontal_pairs = [
        (0, 4), (2, 4),      (3, 7), (5, 7),      (6, 1), (8, 1),     # Cells 0-2
        (9, 13), (11, 13),   (12, 16), (14, 16),  (15, 10), (17, 10), # Cells 3-5
        (18, 22), (20, 22),  (21, 25), (23, 25),  (24, 19), (26, 19), # Cells 6-8
    ]
    # Now, couplings between unit cells parallel to the b-axis (also two per cell)
    vertical_pairs = [
        (1, 9), (2, 9),      (4, 12), (5, 12),    (7, 15), (8, 15),   # Cells 0-2
        (10, 18), (11, 18),  (13, 21), (14, 21),  (16, 24), (17, 24), # Cells 3-5
        (19, 0), (20, 0),    (22, 3), (23, 3),    (25, 6), (26, 6),   # Cells 6-8
    ]

    couplings = []

    rskw = {'dtype':complex, 'order':'F'}
    couplings.extend(Coupling(idx1, idx2, np.eye(3, **rskw), np.array([1., 0, 0])) for (idx1, idx2) in horizontal_pairs)
    couplings.extend(Coupling(idx2, idx1, np.eye(3, **rskw), np.array([-1., 0, 0])) for (idx1, idx2) in horizontal_pairs)
    couplings.extend(Coupling(idx1, idx2, np.eye(3, **rskw), np.array([0., 1, 0])) for (idx1, idx2) in vertical_pairs)
    couplings.extend(Coupling(idx2, idx1, np.eye(3, **rskw), np.array([0., -1, 0])) for (idx1, idx2) in vertical_pairs)
    couplings.extend(Coupling(idx1, idx2, np.eye(3, **rskw), np.array([0., 0, 0])) for (idx1, idx2) in other_pairs)
    couplings.extend(Coupling(idx2, idx1, np.eye(3, **rskw), np.array([0., 0, 0])) for (idx1, idx2) in other_pairs)

    q_mags = 0.5 * np.linspace(0, 1, n_q + 1).reshape(-1, 1)

    # q_vectors = np.concatenate((
    #         q_mags[::-1].reshape(-1, 1) * np.array([1, 0, 1]).reshape(1, -1),
    #         q_mags[1:].reshape(-1, 1) * np.array([0, 0, 1]).reshape(1, -1)
    # ))
    q_vectors = np.concatenate(
        (
            q_mags[::-1].reshape(-1, 1) * np.array([-1, 0, 0]).reshape(1, -1),
            q_mags[1:].reshape(-1, 1) * np.array([1, 1, 0]).reshape(1, -1),
        )
    )

    # q_vectors = q_mags.reshape(-1, 1) * np.array([1, 1, 0]).reshape(1, -1)

    return (rotations, magnitudes, q_vectors, couplings)

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    structure = kagome_supercell()
    q_vectors = structure[2]

    indices = np.arange(201)

    label_indices = [0, 100, 200]
    # label_indices = []

    labels = [str(q_vectors[idx, :]) for idx in label_indices]

    energies = spinwave_calculation(*structure)

    energies = [np.sort(energy.real) for energy in energies]

    positive_energies = [energy[energy > 0] for energy in energies]
    min_energy = min([np.min(energy) for energy in positive_energies])
    translated_energies = [energy - min_energy for energy in positive_energies]

    # Note: we get complex data types with real part zero

    plt.plot(indices, translated_energies)
    # plt.plot(indices, [method.value for method in result.method])
    plt.xticks(label_indices, labels)

    #plt.savefig("fig.png")
    plt.show()
