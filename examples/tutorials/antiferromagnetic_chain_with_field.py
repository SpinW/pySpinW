""" Antiferromagnetic chain example with applied magnetic field

Same as antiferromagnetic_chain.py, bit with field """

from pyspinw import *

unit_cell = UnitCell(4,6,6)

sites = generate_structure(unit_cell, positions=[[0,0,0], [0.5,0,0]], spins=[[0,0,1], [0,0,-1]], names=['X', 'Y'])

exchanges = generate_exchanges(sites=sites,
                               bond=1,
                               exchange_type=HeisenbergExchange,
                               j=1,
                               direction_filter=filter([1,0,0], symmetric=True))

anisotropies = axis_anisotropies(sites, -0.1)

hamiltonian = Hamiltonian(sites, exchanges, anisotropies)

hamiltonian.print_summary()

path = Path([[0,0,0], [2,0,0]])
hamiltonian.spaghetti_plot(path, field=[0,0,7])
