""" Square-lattice Antiferromagnet example

Reproduces Tutorial 11: https://spinw.org/tutorials/11tutorial"""

from pyspinw import *

import matplotlib.pyplot as plt

unit_cell = UnitCell(3, 3, 9)

sites = [LatticeSite(0, 0, 0, 1, 0, 0, name="X")]
k = CommensuratePropagationVector(0.5, 0.5, 0)
s = Structure(sites, unit_cell, supercell=SummationSupercell(propagation_vectors=[k]))

j1 = generate_exchanges(sites, unit_cell, bond=1, exchange_type=HeisenbergExchange, j=59.65)
j2 = generate_exchanges(sites, unit_cell, bond=2, exchange_type=HeisenbergExchange, j=-3.75)
j3 = generate_exchanges(sites, unit_cell, bond=3, exchange_type=HeisenbergExchange, j=1)
exchanges = j1 + j2 + j3

hamiltonian = Hamiltonian(s, exchanges)
hamiltonian.print_summary()

path = Path([[3/4, 1/4, 0], [1/2, 1/2, 0], [1/2, 0, 0], [3/4, 1/4, 0], [1, 0, 0], [1/2, 0, 0]], n_points_per_segment=51)

fig = hamiltonian.spaghetti_plot(path, show=False)
fig.axes[0].set_ylim(0, 300)
plt.show()
