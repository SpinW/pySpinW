## title: Yttrium Iron Garnet Example

# This tutorial uses a CIF file for specifying the crystal structure, [download it here](YIG_10K.cif)
## file: YIG_10K.cif

# Import numpy, pyspinw and load in the structure
import numpy as np
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

"""Neutron scattering [9] indicates
that nearest-neighbor exchange dominates and we adopt
Jad ¼ −9.60 × 10−21 J, Jdd ¼ −3.24 × 10−21 J, and Jaa ¼
−0.92 × 10−21 J [3], where the subscripts refer to the
majority (d) and minority (a) spins."""

Jad = 9.60e-21 * J_to_meV
Jdd = 3.24e-21 * J_to_meV
Jaa = 0.92e-21 * J_to_meV

# Find just a single ion in for each of the "a" and "d" type sublattices
base_site_a = structure.site_by_name("Fe1")
base_site_d = structure.site_by_name("Fe2")

# Fe1 to Fe2
neighbour, offset = structure.one_neighbour(base_site_a, name_filter="Fe2")
exchange_ad = HeisenbergExchange(base_site_a, neighbour, cell_offset=offset, j=Jad, name="Jad", color=(.6,.1,.6))

# Fe1 to Fe1
neighbour, offset = structure.one_neighbour(base_site_a, name_filter="Fe1")
exchange_aa = HeisenbergExchange(base_site_a, neighbour, cell_offset=offset, j=Jaa, name="Jaa", color=(.9,.1,.1))

# Fe2 to Fe2
neighbour, offset = structure.one_neighbour(base_site_d, name_filter="Fe2")
exchange_dd = HeisenbergExchange(base_site_d, neighbour, cell_offset=offset, j=Jdd, name="Jdd", color=(.1,.1,.9))

exchanges = [exchange_ad, exchange_aa, exchange_dd]

hamiltonian = Hamiltonian(structure, exchanges).symmetry_filled()

field = np.array([0, 0, 0.01])

view(hamiltonian)

