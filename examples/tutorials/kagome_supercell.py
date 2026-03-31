""" Kagome 3x3 Antiferromagnet example

Reproduces Tutorial 8: https://spinw.org/tutorials/08tutorial"""

from pyspinw import *
from pyspinw.legacy.genmagstr import genmagstr

unit_cell = UnitCell(6, 6, 40, gamma=120)

x = LatticeSite(0.5, 0,   0, 0, 1, 0, name="X")
y = LatticeSite(0,   0.5, 0, 0, 1, 0, name="Y")
z = LatticeSite(0.5, 0.5, 0, -1, -1, 0, name="Z")
s = genmagstr([x, y, z], unit_cell, magnitude=[1,1,1], mode='helical', k=[-1./3, -1./3, 0], n=[0, 0, 1], unit='lu')

exchanges = generate_exchanges(sites=[x, y, z],
                               unit_cell=unit_cell,
                               max_distance=3.1,
                               coupling_type=HeisenbergCoupling,
                               j=1)

hamiltonian = Hamiltonian(s, exchanges)

hamiltonian.print_summary()


path = Path([[-0.5,0,0], [0,0,0], [0.5,0.5,0]])

import matplotlib.pyplot as plt
fig = hamiltonian.spaghetti_plot(path, show=False)
fig.axes[0].set_ylim(0, 3)
fig.axes[1].set_ylim(0, 1)
plt.show()
