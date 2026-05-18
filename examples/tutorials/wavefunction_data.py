""" Kagome 3x3 Antiferromagnet example

Reproduces Tutorial 8: https://spinw.org/tutorials/08tutorial"""
import numpy as np
import matplotlib.pyplot as plt

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

for i in range(3,4):
    n1 = i
    n2 = i+1

    w1 = wavefunctions[n1]
    w2 = wavefunctions[n2]

    mag1 = np.abs(w1.conj().T @ w1)
    mag2 = np.abs(w2.conj().T @ w2)

    mag_cross = np.abs(w1.conj().T @ w2)

    plt.figure()
    #
    # plt.subplot(2,2,1)
    # plt.pcolor(mag1/scaling_matrix(mag1, mag1))
    # plt.title(f"Correlation matrix {n1}")
    #
    # plt.subplot(2,2,2)
    # plt.pcolor(mag2/scaling_matrix(mag2, mag2))
    # plt.title(f"Correlation matrix {n2}")
    #
    # plt.subplot(2,2,3)

    connection_matrix = mag_cross/scaling_matrix(mag2, mag1)
    above_threshold = connection_matrix > 0.05


    plt.pcolor(above_threshold)
    plt.title(f"Cross correlation matrix {n1} - {n2}")
    plt.colorbar()


plt.show()
