""" Problem set 2 - Supercells """

# Part 1: Make a single lattice unit antiferromagnetic chain with a single unit cell
#         You can use the TransformationSupercell, SummationSupercell, RotationSupercell to do this
#         There are helper functions, e.g. generate_helical_structure
#         Bonus points for all three
#
# Part 2: Make a square (2D) antiferromagnetic lattice with the following Heisenberg exchanges:
#         59.65 between nearest neighbours
#         -3.75 between second-nearest neighbours
#          1.00 between third-nearest neighbours
#
#         How important is the 3rd nearest neighbour interaction?


from pyspinw import *

unit_cell = UnitCell(3, 8, 8)

sites = [LatticeSite(0, 0, 0, 0, 1, 0, name="Bob")]

#
#  FILL IN HERE
#


hamiltonian = Hamiltonian(sites, exchanges)
hamiltonian.print_summary()

path = Path([[0,0,0], [1,0,0]])
hamiltonian.spaghetti_plot(path, scale='log')
