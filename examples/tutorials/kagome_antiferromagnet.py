""" Kagome Antiferromagnet example

Reproduces Tutorial 7: https://spinw.org/tutorials/07tutorial"""

from pyspinw import *
from pyspinw.legacy.genmagstr import genmagstr

unit_cell = UnitCell(6, 6, 10, gamma=120)

s = generate_structure(unit_cell, positions=[[0.5,0,0], [0,0.5,0], [0.5,0.5,0]], moments=[[1,2,0], [-2,-1,0], [1,-1,0]],
                       names=['X','Y','Z'], magnitudes=[1,1,1], moments_unit='lu')

j1 = generate_exchanges(sites=s, bond=1, coupling_type=HeisenbergCoupling, j=1)
j2 = generate_exchanges(sites=s, bond=2, coupling_type=HeisenbergCoupling, j=0.11)
exchanges = j1 + j2

hamiltonian = Hamiltonian(s, exchanges)

hamiltonian.print_summary()

path = Path([[-0.5,0,0], [0,0,0], [0.5,0.5,0]])
hamiltonian.spaghetti_plot(path)
