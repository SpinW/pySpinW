""" Problem set 1 - Simple chains """

from pyspinw import *

unit_cell = UnitCell(1,1,1)

# define a spin in z direction
only_site = LatticeSite(0, 0, 0, 0, 0, 1 , name="X")

s = Structure([only_site], unit_cell=unit_cell)

exchanges = ... # PROBLEM 1: Use HeisebergExchange or generate_exchanges to specify spins
                # HINT: think about the path

hamiltonian = Hamiltonian(s, exchanges)

path = Path([[0,0,0], [1,0,0]])

hamiltonian.spaghetti_plot(path, dE=0.4)

# PROBLEM 2: Create an antiferromagnetic chain with two sites per unit cell
