""" Kagome 3x3 Antiferromagnet example

Reproduces Tutorial 8: https://spinw.org/tutorials/08tutorial"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix

from scipy.sparse.csgraph import min_weight_full_bipartite_matching

from pyspinw import *

unit_cell = UnitCell(3, 3, 4, gamma=120)

sites = generate_helical_structure(unit_cell, positions=[[0,0,0]], spins=[[0, 1, 0]], magnitudes=[3. / 2], names=['X'],
                                   perpendicular=[0,0,1], propagation_vector=[1./3., 1./3., 0])

exchanges = generate_exchanges(sites=sites,
                               max_distance=3.1,
                               exchange_type=HeisenbergExchange,
                               j=1)

anisotropies = axis_anisotropies(sites, 0.2)
hamiltonian = Hamiltonian(sites, exchanges, anisotropies)

hamiltonian.print_summary()

path = Path([[-0.5,0,0], [0,0,0], [0.5,0.5,0]])

energy, intensity, _, wavefunctions = hamiltonian.spinwave_calculation(path, save_wavefunctions=True, use_rust=False, use_rotating=False)

wavefunctions = [wavefunction for wavefunction in wavefunctions]

def scaling_matrix(a, b):
    return np.sqrt(np.diag(a)[None, :] * np.diag(b)[:, None])

def find_ordering(wavefunctions):
    wavefunctions = list(wavefunctions)

    order = np.arange(0, wavefunctions[0].shape[0], dtype=int)
    orderings = [order]

    for w1, w2 in zip(wavefunctions, wavefunctions[1:]):

        mag1 = np.abs(w1.conj().T @ w1)
        mag2 = np.abs(w2.conj().T @ w2)

        mag_cross = np.abs(w1.conj().T @ w2)

        connection_matrix = mag_cross / scaling_matrix(mag2, mag1)

        below_threshold = connection_matrix < 0.05
        connection_matrix[below_threshold] = 0.0
        graph = coo_matrix(connection_matrix)

        # to_inds, from_inds = min_weight_full_bipartite_matching(graph, maximize=True)
        from_inds, to_inds = min_weight_full_bipartite_matching(graph, maximize=True)
        new_order = np.empty_like(order, dtype=int)
        new_order[to_inds] = order[from_inds]

        order = new_order
        orderings.append(order)

    return np.array(orderings, dtype=int)

def apply_orderings(data, orderings):
    data = np.array(data)

    # dimension zero is each q point entry

    for i in range(data.shape[0]):
        data[i, :] = data[i, orderings[i, :]]

    return data


orderings = find_ordering(wavefunctions)
ordered = apply_orderings(energy, orderings)

print(ordered.shape)


plt.figure()
for i in range(ordered.shape[1]):
    plt.plot(np.array(energy)[:, i])


plt.figure()
for i in range(ordered.shape[1]):
    plt.plot(ordered[:, i])



# plt.show()
n_modes = wavefunctions[0].shape[0]
boson_identity = np.zeros((n_modes, n_modes))
boson_identity[:, n_modes//2:] *= -1

for i in range(3,4):
    n1 = i
    n2 = i+1

    w1 = wavefunctions[n1]
    w2 = wavefunctions[n2]

    mag1 = np.abs(w1.conj() @ boson_identity @ w1.T)
    mag2 = np.abs(w2.conj() @ boson_identity @ w2.T)

    mag_cross = np.abs(w1.conj() @ boson_identity @ w2.T)

    plt.figure()

    plt.subplot(2,2,1)
    plt.pcolor(mag1)
    plt.title(f"Correlation matrix {n1}")

    plt.subplot(2,2,2)
    plt.pcolor(mag2)
    plt.title(f"Correlation matrix {n2}")



    plt.subplot(2,2,3)
    plt.pcolor(mag_cross)

    plt.subplot(2,2,4)

    connection_matrix = mag_cross/scaling_matrix(mag2, mag1)
    above_threshold = connection_matrix > 0.05


    plt.pcolor(connection_matrix)
    plt.title(f"Cross correlation matrix {n1} - {n2}")
    plt.colorbar()


plt.show()
