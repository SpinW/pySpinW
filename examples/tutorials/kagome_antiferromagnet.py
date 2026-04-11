""" Kagome Antiferromagnet example

Reproduces Tutorial 7: https://spinw.org/tutorials/07tutorial"""

from pyspinw import *

unit_cell = UnitCell(6, 6, 10, gamma=120)

s = generate_structure(unit_cell, positions=[[0.5,0,0], [0,0.5,0], [0.5,0.5,0]], spins=[[1,2,0], [-2,-1,0], [1,-1,0]],
                       names=['X','Y','Z'], magnitudes=[1,1,1], spins_unit='lu')

j1 = generate_exchanges(sites=s, bond=1, exchange_type=HeisenbergExchange, j=1)
j2 = generate_exchanges(sites=s, bond=2, exchange_type=HeisenbergExchange, j=0.11)
exchanges = j1 + j2

hamiltonian = Hamiltonian(s, exchanges)

hamiltonian.print_summary()

path = Path([[-0.5,0,0], [0,0,0], [0.5,0.5,0]])
hamiltonian.spaghetti_plot(path)
