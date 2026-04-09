""" Antiferromagnetic chain example

Reproduces Tutorial 2: https://spinw.org/tutorials/02tutorial"""

from pyspinw import *

unit_cell = UnitCell(3, 8, 8)

sites = [LatticeSite(0, 0, 0, 0, 1, 0, name="MCu1")]

sites = generate_helical_structure(unit_cell, positions=[[0,0,0]], spins=[[0, 1, 0]],
                                   perpendicular=[0,0,1], propagation_vector=[0.5, 0, 0], names=["MCu1"])

exchanges = generate_exchanges(sites=sites,
                               max_distance=3.1,
                               exchange_type=HeisenbergExchange,
                               j=1)


hamiltonian = Hamiltonian(s, exchanges)
hamiltonian.print_summary()

path = Path([[0,0,0], [1,0,0]])
hamiltonian.spaghetti_plot(path, scale='log')
