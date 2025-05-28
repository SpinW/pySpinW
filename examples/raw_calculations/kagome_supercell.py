"""Example of a Kagome supercell.

This is the magnet in MATLAB SpinW tutorial 8:
    https://spinw.org/tutorials/08tutorial
"""

import numpy as np

from pyspinw.calculations.spinwave import spinwave_calculation, Coupling


# define our rotation matrices
def rotation(theta):
    """Return the rotation matrix for an angle `theta` in the x-z plane."""
    return -1 * np.array(
        [
            [np.cos(theta), -np.sin(theta), 0],
            [np.sin(theta), np.cos(theta), 0],
            [0, 0, 1],
        ]
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
    def down():
        # ↙️
        return rotation(np.pi / 6)

    def up():
        # ↖
        return rotation(7 * np.pi / 6)

    def right():
        # →
        return rotation(np.pi / 2)

    rotations = np.array(
        [
            down(),
            right(),
            right(),
            up(),
            down(),
            down(),
            right(),
            up(),
            up(),
            up(),
            down(),
            down(),
            right(),
            up(),
            up(),
            down(),
            right(),
            right(),
            right(),
            up(),
            up(),
            down(),
            right(),
            right(),
            up(),
            down(),
            down(),
        ]
    )
    magnitudes = np.array([1] * 27)  # spin 1

    k = 0.5

    horizontal_pairs = [  # couplings which are parallel to the a-axis
        (0, 2), (2, 3), (3, 5), (5, 6), (6, 8),  # bottom horizontal line
        (9, 11), (11, 12), (12, 14), (14, 15), (15, 17),  # middle horizontal line
        (18, 20), (20, 21), (21, 23), (23, 24), (24, 26),  # top horizontal line
        (8, 0), (17, 9), (26, 18),  # 'over the edges' from left
    ]

    vertical_pairs = [  # pairs which are parallel to the b-axis
        (0, 1), (1, 9), (9, 10), (10, 18), (18, 19),  # leftmost vertical line
        (3, 4), (4, 12), (12, 13), (13, 21), (21, 22),  # middle vertical line
        (6, 7), (7, 15), (15, 16), (16, 24), (24, 25),  # rightmost vertical line
        (19, 0), (22, 3), (25, 6),  # 'over the edges' from bottom
    ]

    other_pairs = [  # pairs going diagonally up and to the right
        (1, 11), (11, 13), (13, 23), (23, 25), (8, 1),  # line extending through site 1 (including over edge to 8)
        (2, 4), (4, 14), (14, 16), (16, 26), (19, 2),  # line extending through site 2 (including over edge to 19)
        (5, 7), (7, 17), (22, 5),  # line extending through site 5 (including over edge to 22)
        (10, 20), (20, 22), (17, 10),  # line extending through site 10 (including over edge to 17)
        (25, 8), (26, 19),  # over edges which aren't otherwise part of a line
    ]

    couplings = []

    k = 0.5
    couplings.extend(Coupling(idx1, idx2, rotations[idx1] @ rotations[idx2].T, k*np.array([1, 0, 0])) for (idx1, idx2) in horizontal_pairs)
    couplings.extend(Coupling(idx2, idx1, rotations[idx2] @ rotations[idx1].T, k*np.array([-1, 0, 0])) for (idx1, idx2) in horizontal_pairs)
    couplings.extend(Coupling(idx1, idx2, rotations[idx1] @ rotations[idx2].T, k*np.array([0, 1, 0])) for (idx1, idx2) in vertical_pairs)
    couplings.extend(Coupling(idx2, idx1, rotations[idx2] @ rotations[idx1].T, k*np.array([0, -1, 0])) for (idx1, idx2) in vertical_pairs)
    couplings.extend(Coupling(idx1, idx2, rotations[idx1] @ rotations[idx2].T, k*np.array([1, 1, 0])) for (idx1, idx2) in other_pairs)
    couplings.extend(Coupling(idx2, idx1, rotations[idx2] @ rotations[idx1].T, k*np.array([-1, -1, 0])) for (idx1, idx2) in other_pairs)

    q_mags = 0.5 * np.linspace(0, 1, n_q + 1).reshape(-1, 1)

    # q_vectors = np.concatenate((
    #         q_mags[::-1].reshape(-1, 1) * np.array([1, 0, 1]).reshape(1, -1),
    #         q_mags[1:].reshape(-1, 1) * np.array([0, 0, 1]).reshape(1, -1)
    # ))
    q_vectors = np.concatenate(
        (
            q_mags[::-1].reshape(-1, 1) * np.array([1, 1, 0]).reshape(1, -1),
            q_mags[1:].reshape(-1, 1) * np.array([0, 1, 0]).reshape(1, -1),
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

    energies = spinwave_calculation(*structure).raw_energies

    energies = [np.sort(energy.real) for energy in energies]

    positive_energies = [energy[energy > 0] for energy in energies]
    min_energy = min([np.min(energy) for energy in positive_energies])
    translated_energies = [energy - min_energy for energy in positive_energies]

    # Note: we get complex data types with real part zero

    plt.plot(indices, translated_energies)
    # plt.plot(indices, [method.value for method in result.method])
    plt.xticks(label_indices, labels)

    plt.savefig("fig.png")
