import numpy as np

## title: Yttrium Iron Garnet Example

# This tutorial uses a CIF file for specifying the crystal structure, [download it here](YIG_10K.cif)
## file: YIG_10K.cif

# Import pyspinw and load in the structure

from pyspinw import *

structure = load("YIG_10K.cif")

# This structure has two iron sublattices
# Colour the sublattices associated with Fe1 and Fe2 atoms

for site in structure.sites_by_name("Fe1"):
    site.metadata.color = 1,0,0
    site.spin = 3/2, 0, 0

for site in structure.sites_by_name("Fe2"):
    site.metadata.color = 0,0,1
    site.spin = 3/2, 0, 0

J_to_meV = 6.241509e21

Jad = 9.60e-21 * J_to_meV
Jdd = 3.24e-21 * J_to_meV
Jaa = 0.92e-21 * J_to_meV

exchanges = []
base_site = structure.site_by_name("Fe1")
for i, (energy, name) in enumerate(zip([Jad, Jdd, Jaa], ["Jad", "Jdd", "Jaa"])):
    other_site, offset = structure.one_neighbour(base_site, n=i+1, element="Fe")

    exchanges.append(
        HeisenbergExchange(base_site, other_site, cell_offset=offset, name=name, j=energy))

hamiltonian = Hamiltonian(structure, exchanges).symmetry_filled()

field = np.array([0, 0, 0.01])

view(hamiltonian)

