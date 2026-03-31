""" Kagome Ferromagnet example

Reproduces Tutorial 5: https://spinw.org/tutorials/05tutorial"""

from pyspinw import *

unit_cell = UnitCell(6,6,5, gamma=120)

x = LatticeSite(0.5, 0, 0, 0, 1, 0, name="X")
y = LatticeSite(0, 0.5, 0, 0, 1, 0, name="Y")
z = LatticeSite(0.5, 0.5, 0, 0, 1, 0, name="Z")

sites = [x, y, z]

s = Structure(sites, unit_cell=unit_cell, supercell=TrivialSupercell(scaling=(1,1,1)))


exchanges = generate_exchanges(sites=[x, y, z],
                               unit_cell=unit_cell,
                               max_distance=4.,
                               coupling_type=HeisenbergCoupling,
                               j=-1)

hamiltonian = Hamiltonian(s, exchanges)

hamiltonian.print_summary()

path = Path([[-0.5,0,0], [0,0,0], [0.5,0.5,0]])
hamiltonian.energy_plot(path, show=False)

sample = Powder(hamiltonian)
path1D = Path1D(0.1, 2.5, n_points=100)
sample.show_spectrum(path1D, n_energy_bins=250)
