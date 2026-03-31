""" Antiferromagnetic chain example with applied magnetic field

Same as antiferromagnetic_chain.py, bit with field """

from pyspinw import *

unit_cell = UnitCell(4,6,6)

sites = [LatticeSite(0, 0, 0, 0,0,1, name="X"),
         LatticeSite(0.5,0,0, 0,0, -1, name="Y")]

s = Structure(sites, unit_cell=unit_cell)

exchanges = generate_exchanges(sites=sites,
                               unit_cell=unit_cell,
                               max_distance=2.1,
                               coupling_type=HeisenbergCoupling,
                               j=1,
                               direction_filter=filter([1,0,0], symmetric=True))

anisotropies = axis_anisotropies(sites, -0.1)

hamiltonian = Hamiltonian(s, exchanges, anisotropies)

hamiltonian.print_summary()

path = Path([[0,0,0], [2,0,0]])
hamiltonian.spaghetti_plot(path, field=[0,0,7])
