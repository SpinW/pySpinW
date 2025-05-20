"""Example of a Kagome supercell.

This is the magnet in MATLAB SpinW tutorial 8:
    https://spinw.org/tutorials/08tutorial
"""

import numpy as np

from pyspinw.calculations.spinwave import spinwave_calculation, Coupling

def sym_coupling(idx1, idx2, rot1, rot2, inter_site_vector):
    """Create a symmetric coupling."""
    # the matrix for a coupling is
    # R_{idx1} @ R_{idx2}^{-1}
    # and rotation matrices are orthogonal so this is equal to
    # R_{idx1} @ R_{idx2}^T
    return (Coupling(idx1, idx2, rot1 @ rot2.T, inter_site_vector),
            Coupling(idx2, idx1, rot2 @ rot1.T, -inter_site_vector))

def unit_cell_couplings(cell_number: int, rotations):
    """Create all the couplings within a unit cell."""
    k = 0.5
    c = cell_number*3
    return      [
                 *sym_coupling(c, c+1, rotations[c], rotations[c+1], k*np.array([ 1,  0, 0])),
                 *sym_coupling(c, c+1, rotations[c], rotations[c+1], k*np.array([-1,  0, 0])),
                 *sym_coupling(c, c+2, rotations[c], rotations[c+2], k*np.array([ 1,  1, 0])),
                 *sym_coupling(c, c+2, rotations[c], rotations[c+2], k*np.array([-1, -1, 0])),
                 *sym_coupling(c+1, c+2, rotations[c+1], rotations[c+2], k*np.array([ 0,  1, 0])),
                 *sym_coupling(c+1, c+2, rotations[c+1], rotations[c+2], k*np.array([ 0, -1, 0])),
                 ]

# define our rotation matrices
def rotation(theta):
    """Return the rotation matrix for an angle `theta` in the x-z plane."""
    return np.array([[np.cos(theta), np.sin(theta), 0],
                    [0, 0, -1],
                    [np.sin(theta), -np.cos(theta), 0]])

def kagome_supercell():
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
        return rotation(7*np.pi/6)
    def up():
        # ↖
        return rotation(np.pi/6)
    def right():
        # →
        return rotation(np.pi/2)

    rotations = np.array([down(),
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
                          down()
                          ])
    magnitudes = np.array([1]*27)  # spin 1

    k = 0.5

    couplings = []

    for cell in range(0, 8):
        couplings.extend(unit_cell_couplings(cell, rotations))

    n_q = 101
    q_mags = 0.5*np.linspace(0, 1, n_q).reshape(-1, 1)

    # q_vectors = np.concatenate((
    #         q_mags[::-1].reshape(-1, 1) * np.array([1, 0, 1]).reshape(1, -1),
    #         q_mags[1:].reshape(-1, 1) * np.array([0, 0, 1]).reshape(1, -1)
    # ))
    q_vectors = np.concatenate((
            q_mags[::-1].reshape(-1, 1) * np.array([1, 1, 0]).reshape(1, -1),
            q_mags[1:].reshape(-1, 1) * np.array([0, 1, 0]).reshape(1, -1)
    ))

    # q_vectors = q_mags.reshape(-1, 1) * np.array([1, 1, 0]).reshape(1, -1)

    indices = np.arange(201)

    label_indices = [0, 100, 200]
    # label_indices = []

    labels = [str(q_vectors[idx,:]) for idx in label_indices]

    structure = (rotations, magnitudes, q_vectors, couplings)

    return structure, indices, labels, label_indices

if __name__ == "__main__":

    import matplotlib.pyplot as plt

    structure, indices, labels, label_indices = kagome_supercell()
    energies = spinwave_calculation(*structure).raw_energies

    energies = [np.sort(energy.real) for energy in energies]

    positive_energies = [energy[energy>0] for energy in energies]
    min_energy = min([np.min(energy) for energy in positive_energies])
    translated_energies = [energy - min_energy for energy in positive_energies]


    # Note: we get complex data types with real part zero

    plt.plot(indices, translated_energies)
    # plt.plot(indices, [method.value for method in result.method])
    plt.xticks(label_indices, labels)

    plt.savefig("fig.png")
