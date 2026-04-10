""" Spin wave spectrum of Ba3NbFe3Si2O14 """
from pyspinw import *

unit_cell = UnitCell(8.539,8.539, 5.2414, 90, 90, 120)
sites = [LatticeSite(0, 0.24964, 1/2, 0, 0, 5/2, name="Fe3")]
sg = spacegroup("P 3 2 1")


structure = Structure(sites, unit_cell, sg, supercell=TiledSupercell(scaling=(2,2,1)))
# structure = Structure(sites, unit_cell, sg)

exchanges = generate_exchanges(structure.sites, unit_cell=unit_cell, max_distance=5) + \
            generate_exchanges(structure.sites, unit_cell=unit_cell, min_distance=5, max_distance=8,
                               direction_filter=filter([0,0,1], True))

hamiltonian = Hamiltonian(structure, exchanges)

hamiltonian.print_summary()


view(hamiltonian)

#exchanges = generate_exchanges(sites)