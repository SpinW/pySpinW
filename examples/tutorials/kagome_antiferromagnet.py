""" Kagome Antiferromagnet example

Reproduces Tutorial 7: https://spinw.org/tutorials/07tutorial"""

from pyspinw import *
from pyspinw.legacy.genmagstr import genmagstr

unit_cell = UnitCell(6, 6, 10, gamma=120)

x = LatticeSite(0.5, 0, 0, 1, 2, 0, name="X")
y = LatticeSite(0, 0.5, 0, -2, -1, 0, name="Y")
z = LatticeSite(0.5, 0.5, 0, 1, -1, 0, name="Z")
s = genmagstr([x, y, z], unit_cell, magnitude=[1,1,1], mode='tile', unit='lu')

j1 = generate_exchanges(sites=[x, y, z], unit_cell=unit_cell, min_distance=0, max_distance=4.1, coupling_type=HeisenbergCoupling, j=1)
j2 = generate_exchanges(sites=[x, y, z], unit_cell=unit_cell, min_distance=5, max_distance=5.2, coupling_type=HeisenbergCoupling, j=0.11)
exchanges = j1 + j2

hamiltonian = Hamiltonian(s, exchanges)

hamiltonian.print_summary()

path = Path([[-0.5,0,0], [0,0,0], [0.5,0.5,0]])
hamiltonian.spaghetti_plot(path)
