""" Antiferromagnetic chain example

Reproduces Tutorial 2: https://spinw.org/tutorials/02tutorial"""

from pyspinw import *

unit_cell = UnitCell(3, 8, 8)

sites = [LatticeSite(0, 0, 0, 0, 1, 0, name="MCu1")]

k = CommensuratePropagationVector(0.5, 0, 0)
s = Structure(sites, unit_cell=unit_cell, supercell=SummationSupercell(propagation_vectors=[k]))

exchanges = generate_exchanges(sites=sites,
                               unit_cell=unit_cell,
                               max_distance=3.1,
                               coupling_type=HeisenbergCoupling,
                               j=1)


hamiltonian = Hamiltonian(s, exchanges)
hamiltonian.print_summary()

path = Path([[0,0,0], [1,0,0]])
hamiltonian.spaghetti_plot(path, scale='log')
