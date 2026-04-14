""" Ferromagnetic chain example """
import numpy as np
import matplotlib.pyplot as plt

from pyspinw import *

unit_cell = UnitCell(1,1,1)

only_site = LatticeSite(0, 0, 0, 0,0,1, name="X")

s = Structure([only_site], unit_cell=unit_cell)

exchanges = generate_exchanges(sites=[only_site],
                               unit_cell=unit_cell,
                               max_distance=1.1,
                               exchange_type=HeisenbergExchange,
                               j=-1,
                               direction_filter=filter([1,0,0]))

hamiltonian = Hamiltonian(s, exchanges)

path_x = Path([[0,0,0], [1,0,0]])
path_y = Path([[0,0,0], [0,1,0]])
path_z = Path([[0,0,0], [0,0,1]])

path_xy = Path([[0,0,0], [1,1,0]])
path_xz = Path([[0,0,0], [1,0,1]])

# print("Hamiltonian Spaghetti")
#
# hamiltonian.spaghetti_plot(path_x)
# hamiltonian.spaghetti_plot(path_y)
# hamiltonian.spaghetti_plot(path_z)
#
# hamiltonian.spaghetti_plot(path, dE=0.4)

#
# print("Single Crystal Spaghetti")
#
# plt.figure("Single Crystal")
# single_crystal = SingleCrystal(hamiltonian)
# plt.subplot(2,2,1)
# single_crystal.spaghetti_plot(path_x, dE=0.4, show=False, new_figure=False)
# plt.subplot(2,2,2)
# single_crystal.spaghetti_plot(path_y, dE=0.4, show=False, new_figure=False)
# plt.subplot(2,2,3)
# single_crystal.spaghetti_plot(path_z, dE=0.4, show=False, new_figure=False)
# plt.subplot(2,2,4)
# single_crystal.spaghetti_plot(path, dE=0.4, show=False, new_figure=False)

print("Twin Spaghetti")

twin = Twin(hamiltonian, np.array([[0,1,0], [-1,0,0], [0,0,1]], dtype=float), 0.2)




for p in [path_x, path_y, path_z, path_xy, path_xz]:
    plt.figure(f"Twin {p}")
    twin.spaghetti_plot(p, dE=0.4, show=False, new_figure=False)

plt.show()

twin.spaghetti_plot(p)