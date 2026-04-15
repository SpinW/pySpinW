from pyspinw import *

unit_cell = UnitCell(1,1,1)

# Define a spin in z direction
only_site = LatticeSite(0, 0, 0, 0, 0, 1 , name="X")

# Create a structure
s = Structure([only_site], unit_cell=unit_cell)

exchanges = [
    HeisenbergExchange(only_site, only_site, cell_offset=(1, 0, 0), j=-1)]

hamiltonian = Hamiltonian(s, exchanges)

path = Path([[0,0,0], [1,0,0]])

hamiltonian.spaghetti_plot(path, dE=0.4)

# PROBLEM 2: Create an antiferromagnetic chain with two sites per unit cell

# One possible solution

x = LatticeSite(0,0,0,0,0, 1)
y = LatticeSite(0.5,0,0,0,0,-1)

s = Structure([x, y], unit_cell=unit_cell)

exchanges = [
    HeisenbergExchange(x, y, cell_offset=(0, 0, 0), j=1),
    HeisenbergExchange(y, x, cell_offset=(1, 0, 0), j=1)]

hamiltonian = Hamiltonian(s, exchanges)

hamiltonian.spaghetti_plot(path)
