""" Antiferromagnetic chain example

Reproduces Tutorial 2: https://spinw.org/tutorials/02tutorial
"""
## title: Antiferromagnetic chain with rotating-frame calculation
## reproduces: 2

# Import pyspinw
from pyspinw import *

# Create a 3x8x8 unit cell
unit_cell = UnitCell(3, 8, 8)

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

# You can get a summary of the Hamiltonian using `print_summary` (also works on `Structure`)
## capture-stdout
hamiltonian.print_summary()
## capture-end

# View it using view
## skip
view(hamiltonian)

## image: structure.png
## snapshot(hamiltonian, filename="structure.png")

# Now we can plot a spaghetti diagram, first we generate a path though the lattice
# from $q=0$ to $1$ reciprocal lattice unit in $x$.
path = Path([[0,0,0], [1,0,0]])

# Show a plot
## skip: 1
hamiltonian.spaghetti_plot(path, scale='log')

## image: spaghetti_plot.png
## fig = hamiltonian.spaghetti_plot(path, scale='log', show=False)
## fig.savefig("spaghetti_plot.png")