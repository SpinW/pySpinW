""" Ferromagnetic chain example """

from pyspinw import *


unit_cell = UnitCell(1,1,1)

only_site = LatticeSite(0, 0, 0, 0,0,1, name="X")

s = Structure([only_site], unit_cell=unit_cell)

exchanges = generate_exchanges(sites=[only_site],
                               unit_cell=unit_cell,
                               max_distance=1.1,
                               exchange_type=HeisenbergExchange,
                               j=-1,
                               direction_filter=filter([1,0,0]))

hamiltonian = Hamiltonian(s, exchanges)

path = Path([[0,0,0], [1,0,0]])

sample = Powder(hamiltonian)

path1D = Path1D(0.01, 1)
sample.show_spectrum(path1D, n_energy_bins=100, n_samples=500, energy_stddev=0.4)