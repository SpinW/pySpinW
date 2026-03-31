""" Antiferromagnetic chain example """
from pickletools import optimize

import numpy as np
from scipy.optimize import minimize_scalar, Bounds, minimize, differential_evolution, least_squares

from pyspinw.coupling import HeisenbergCoupling
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.path import Path1D
from pyspinw.sample import Powder
from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structures import Structure

unit_cell = UnitCell(3, 8, 8)

x = LatticeSite(0, 0, 0, 0, 1, 0, name="X")
y = LatticeSite(0.5, 0, 0, 0, 1, 0, name="Y")

j1 = HeisenbergCoupling(x, y, j=1, cell_offset=(0,0,0), name="J1")
j2 = HeisenbergCoupling(y, x, j=1, cell_offset=(0,1,0), name="J2")

sites = [x, y]
exchanges = [j1, j2]

s = Structure(sites, unit_cell)

hamiltonian = Hamiltonian(s, exchanges)

sample = Powder(hamiltonian)

path1D = Path1D(0.01, 1, n_points=20)

n_energy = 20
n_samples = 500

sample.show_spectrum(path1D, n_energy_bins=n_energy, n_samples=n_samples, energy_stddev=0.4)

spectrum = sample.parameterized_spectrum(
    parameters=["J.j"],
    path=path1D,
    n_energy_bins=n_energy,
    n_samples=n_samples,
    energy_stddev=0.4,
    find_ground_state_with={"fixed": [x], "verbose": False})

target = spectrum(1.2)

def objective(x):
    return np.sum((target - spectrum(x))**2) / 1e10

def objective_lin(x):
    return (target - spectrum(x)).reshape(-1)

def callback(intermediate_result):
    """ Callback to print out what is going on"""
    print(intermediate_result)

x = np.linspace(0.9, 1.3, 31)
import matplotlib.pyplot as plt
plt.plot(x, [objective(xx) for xx in x])
plt.show()


print("Optimising...")

import time
t0 = time.time()

bounds = Bounds([0.5], [1.5])
solution = least_squares(objective_lin, x0=1.0, bounds=bounds, callback=callback, diff_step=1e-2)

print(solution)
print(time.time() - t0)

