import numpy as np
from pyspinw.calculations.spinwave import spinwave_calculation, Coupling

def kagome_ferromagnet(n_q = 100):

    """
    Basic ferromagnet
    """

    # Three sites, otherwise identical
    rotations = np.array([np.eye(3) for _ in range(3)])
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
                 Coupling(0, 1, np.eye(3), inter_site_vector=k*np.array([ 1,  0, 0])),
                 Coupling(0, 1, np.eye(3), inter_site_vector=k*np.array([-1,  0, 0])),
                 Coupling(1, 0, np.eye(3), inter_site_vector=k*np.array([ 1,  0, 0])),
                 Coupling(1, 0, np.eye(3), inter_site_vector=k*np.array([-1,  0, 0])),
                 Coupling(0, 2, np.eye(3), inter_site_vector=k*np.array([ 1,  1, 0])),
                 Coupling(0, 2, np.eye(3), inter_site_vector=k*np.array([-1, -1, 0])),
                 Coupling(2, 0, np.eye(3), inter_site_vector=k*np.array([ 1,  1, 0])),
                 Coupling(2, 0, np.eye(3), inter_site_vector=k*np.array([-1, -1, 0])),
                 Coupling(1, 2, np.eye(3), inter_site_vector=k*np.array([ 0,  1, 0])),
                 Coupling(1, 2, np.eye(3), inter_site_vector=k*np.array([ 0, -1, 0])),
                 Coupling(2, 1, np.eye(3), inter_site_vector=k*np.array([ 0,  1, 0])),
                 Coupling(2, 1, np.eye(3), inter_site_vector=k*np.array([ 0, -1, 0]))
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

    return (rotations, magnitudes, q_vectors, couplings)

if __name__ == "__main__":

    import matplotlib.pyplot as plt

    structure = kagome_ferromagnet()
    q_vectors = structure[2]

    indices = np.arange(201)

    label_indices = [0, 100, 200]
    # label_indices = []

    labels = [str(q_vectors[idx,:]) for idx in label_indices]

    result = spinwave_calculation(*structure)

    plt.plot(indices, result.raw_energies)
    # plt.plot(indices, [method.value for method in result.method])
    plt.xticks(label_indices, labels)

    # Compare with tutorial 5, 3rd last figure (https://spinw.org/tutorial5_05.png)
    plt.show()
