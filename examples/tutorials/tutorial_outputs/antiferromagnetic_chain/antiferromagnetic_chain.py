""" Antiferromagnetic chain example

Reproduces Tutorial 2: https://spinw.org/tutorials/02tutorial
"""


# In this tutorial we will then take three different approaches to calculating the magnon spectrum
#  - Creating unit cell containing with two spins pointing in different directions.
#  - Using a unit cell with a single spin, performing a calculation that *explicitly* creates two unit cells
#    containing opposing spins.
#  - Using the "rotating frame" optimisation, which *implicitly* defines two cells. This has advantages
#    in terms of speed and clarity for systems where it can be applied.
#
# As in the previous tutorial we start by importing pyspinw
from pyspinw import *


# We create a magnetic structure with a 2x1x1 unit cell and two atoms - at 0.5 and 1.5 angstroms in $x$,
# facing in different directions in $z$.

two_site_unit_cell = UnitCell(2, 1, 1)
two_site_site_1 = LatticeSite(0.25, 0.5, 0.5, sz=1, name="a")
two_site_site_2 = LatticeSite(0.75, 0.5, 0.5, sz=-1, name="b")
two_site_structure = Structure([two_site_site_1, two_site_site_2], two_site_unit_cell)

# View it
view(two_site_structure)

# Like in the previous example, there are different ways we can set up the exchanges, if we do it explicitly we
# have

two_site_exchanges = [
    HeisenbergExchange(two_site_site_1, two_site_site_2, j=1),
    HeisenbergExchange(two_site_site_2, two_site_site_1, cell_offset=(1,0,0), j=1)
]

# Or using `generate_exchanges`, which should give the same result
two_site_exchanges = generate_exchanges([two_site_site_1, two_site_site_2],
                               unit_cell=two_site_unit_cell,
                               max_distance=1.1,
                               direction_filter=filter([1,0,0], symmetric=True))

# If we check to see if it is the same we find that it is not identical to the hand-crafted version
# but it is equivalent, we have exchanges between a and b, one of which is offset by one cell in the
# $x$ (lattice a) direction
for exchange in two_site_exchanges:
    print(exchange)

# Now we build a hamiltonian
two_site_hamiltonian = Hamiltonian(two_site_structure, two_site_exchanges)
# And plot the dispersion, here we specify the path in angstroms to make it comparable with different unit cells.
path = Path([(0, 0, 0), (1, 0, 0)], convert_to_lattice_units_with=two_site_unit_cell)
two_site_hamiltonian.spaghetti_plot(path)

unit_cell = UnitCell(1,1,1)

# The following generates a 2x1x1 supercell in which the spins alternate in y (period 2 rotation around z)
structure = generate_helical_structure(unit_cell, positions=[[0,0,0]], spins=[[0, 1, 0]],
                                   perpendicular=[0,0,1], propagation_vector=[0.5, 0, 0], names=["MCu1"])


# Generate Heisenberg exchanges based on distance, this will only be in x because of the shape of the unit cell
exchanges = generate_exchanges(sites=structure,
                               max_distance=3.1,
                               exchange_type=HeisenbergExchange,
                               j=1)

# Build the Hamiltonian
hamiltonian = Hamiltonian(structure, exchanges)

# And create a spaghetti plot. We will be comparing structures with different unit cells, so,
path = Path([(0,0,0), (1,0,0)])

# You can get a summary of the Hamiltonian using `print_summary` (also works on `Structure`)
hamiltonian.print_summary()

# View it using view
view(hamiltonian)


# Now we can plot a spaghetti diagram, first we generate a path though the lattice
# from $q=0$ to $1$ reciprocal lattice unit in $x$.
path = Path([[0,0,0], [1,0,0]])

# Show a plot
hamiltonian.spaghetti_plot(path, scale='log')

