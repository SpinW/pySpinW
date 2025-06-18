"""Kagome antiferromagnet example.

See https://spinw.org/tutorials/07tutorial.
"""
import numpy as np

try:
    from pyspinw.rust import spinwave_calculation, Coupling
except ModuleNotFoundError:
    from pyspinw.calculations.spinwave import spinwave_calculation, Coupling

# define our rotation matrices
def rotation(theta):
    """Calculates the rotation matrix for a given x-y plane angle theta.

    Returns a matrix whose columns are [e1==real(zed) e2==imag(zed) e3=eta]
    which serves to rotate spins into FM alignment according to
    Toth and Lake JPCM 27 16602 (2016), for a spin in the x-y plane rotated
    at an angle theta to the x-axis.
    """
    """ # From Matlab SpinW
    e3 = np.array([np.sin(theta), np.cos(theta), 0.0])
    e2 = np.cross(e3, np.array([1.0, 0.0, 0.0]))
    if np.linalg.norm(e2) < 1e-6:
        e2 = np.array([0.0, 0.0, 1.0])
    e2 = e2 / np.linalg.norm(e2)
    e1 = np.cross(e2, e3)
    return np.array([e1, e2, e3]).T
    """
    # From Lukas
    return np.array(
        [
            [np.cos(theta), 0, np.sin(theta)],
            [np.sin(theta), 0, -np.cos(theta)],
            [0, 1, 0],
        ],
        dtype=complex, order='F',
    )

def kagome_antiferromagnet(n_q = 100):
    """Kagome anti-ferromagnet like in tutorial 7."""
    # Three sites, otherwise identical
    rotations = [
        rotation(0),
        rotation(2 * np.pi / 3),
        rotation(-2 * np.pi / 3)
    ]
    magnitudes = np.array([1.0]*3)  # spin-1

    # Do the J1 (nearest neighbour) couplings - using table from Matlab Tutorial 7
    # Run the example until the end and then run:
    # >> AFkagome.table('bond',1:2)
    # And use the values in idx1, idx2, and dl (idx1, idx2 indexed from 1 not 0)
    J1mat = np.eye(3, dtype=complex, order='F') * 1.0
    couplings = [
        Coupling(2, 0, J1mat, np.array([0., 1, 0])),
        Coupling(0, 1, J1mat, np.array([0., -1, 0])),
        Coupling(1, 2, J1mat, np.array([0., 0, 0])),
        Coupling(2, 0, J1mat, np.array([0., 0, 0])),
        Coupling(0, 1, J1mat, np.array([1., 0, 0])),
        Coupling(1, 2, J1mat, np.array([-1., 0, 0])),
    ]
    # Also need the conjugate coupling
    couplings.extend([
        Coupling(0, 2, J1mat, np.array([0., -1, 0])),
        Coupling(1, 0, J1mat, np.array([0., 1, 0])),
        Coupling(2, 1, J1mat, np.array([0., 0, 0])),
        Coupling(0, 2, J1mat, np.array([0., 0, 0])),
        Coupling(1, 0, J1mat, np.array([-1., 0, 0])),
        Coupling(2, 1, J1mat, np.array([1., 0, 0])),
    ])

    # Do the J2 (next-nearest-neigbour) couplings
    J2mat = np.eye(3, dtype=complex, order='F') * 0.11
    couplings.extend([
        Coupling(0, 1, J2mat, np.array([1., -1, 0])),
        Coupling(1, 2, J2mat, np.array([0., 1, 0])),
        Coupling(2, 0, J2mat, np.array([-1., 0, 0])),
        Coupling(0, 1, J2mat, np.array([0., 0, 0])),
        Coupling(1, 2, J2mat, np.array([-1., -1, 0])),
        Coupling(2, 0, J2mat, np.array([1., 1, 0])),
    ])
    couplings.extend([
        Coupling(1, 0, J2mat, np.array([-1., 1, 0])),
        Coupling(2, 1, J2mat, np.array([0., -1, 0])),
        Coupling(0, 2, J2mat, np.array([1., 0, 0])),
        Coupling(1, 0, J2mat, np.array([0., 0, 0])),
        Coupling(2, 1, J2mat, np.array([1., 1, 0])),
        Coupling(0, 2, J2mat, np.array([-1., -1, 0])),
    ])

    q_mags = 0.5*np.linspace(0, 1, n_q + 1).reshape(-1, 1)

    q_vectors = np.concatenate((
            q_mags[::-1].reshape(-1, 1) * np.array([-1, 0, 0]).reshape(1, -1),
            q_mags[1:].reshape(-1, 1) * np.array([1, 1, 0]).reshape(1, -1)
    ))

    return (rotations, magnitudes, q_vectors, couplings)

if __name__ == "__main__":

    import matplotlib.pyplot as plt

    structure = kagome_antiferromagnet()
    q_vectors = structure[2]

    indices = np.arange(201)

    label_indices = [0, 100, 200]
    # label_indices = []

    labels = [str(q_vectors[idx,:]) for idx in label_indices]

    energies = spinwave_calculation(*structure)

    # Ignore imaginary energies (we shouldn't get any here...)
    energies = [np.sort(energy.real) for energy in energies]

    positive_energies = [energy[energy>0] for energy in energies]

    plt.plot(indices, positive_energies)
    # plt.plot(indices, [method.value for method in result.method])
    plt.xticks(label_indices, labels)

    # Compare with tutorial 7, 2nd last figure (https://spinw.org/tutorial7_05.png)
    # It looks slightly assymmetric because Matlab-SpinW adjusts the aspect-ratio
    # so that each unit along the x-axis is the same |Q| value whatever the direction
    plt.show()
