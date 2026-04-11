""" Kagome 3x3 Antiferromagnet example

Reproduces Tutorial 8: https://spinw.org/tutorials/08tutorial"""

from pyspinw import *

unit_cell = UnitCell(6, 6, 40, gamma=120)

s = generate_helical_structure(unit_cell, positions=[[0.5,0,0], [0,0.5,0], [0.5,0.5,0]],
                               spins=[[0, 1, 0], [0, 1, 0], [-1, -1, 0]], magnitudes=[1, 1, 1], names=['X', 'Y', 'Z'],
                               spins_unit='lu', perpendicular=[0,0,1], propagation_vector=[-1./3., -1./3., 0])

exchanges = generate_exchanges(sites=s,
                               max_distance=3.1,
                               exchange_type=HeisenbergExchange,
                               j=1)

hamiltonian = Hamiltonian(s, exchanges)

hamiltonian.print_summary()

path = Path([[-0.5,0,0], [0,0,0], [0.5,0.5,0]])

import matplotlib.pyplot as plt
fig = hamiltonian.spaghetti_plot(path, show=False, use_rotating=False)
fig.axes[0].set_ylim(0, 3)
plt.show()
