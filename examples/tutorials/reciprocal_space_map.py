"""Reciprocal space map with energy integration

Demonstrates how to plot a 2D reciprocal space map showing spin wave
intensity integrated over an energy window using the Slice class.
"""

from pyspinw import *

# Create a simple ferromagnetic chain
unit_cell = UnitCell(1, 1, 1)
only_site = LatticeSite(0, 0, 0, 0, 0, 1, name="X")
s = Structure([only_site], unit_cell=unit_cell)

exchanges = generate_exchanges(
    sites=[only_site],
    unit_cell=unit_cell,
    bond=1,
    j=-1
)

hamiltonian = Hamiltonian(s, exchanges)

# Define a 2D slice in reciprocal space
q_slice = Slice(
    origin=[0, 0, 0],
    axis_1=[1, 0, 0],
    axis_2=[0, 1, 0],
    n_a=50,
    n_b=50,
    e_min=0,
    e_max=2,
    padding=0.2
)

hamiltonian.intensity_map(q_slice)
