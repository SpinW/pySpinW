from pyspinw import *

sg = spacegroup("p4/mmm")

x=LatticeSite(0.5, 0.5, 0.5, name="X")
y=LatticeSite(0.0, 0.5, 0.5, name="Y")


cell = UnitCell(1,1,1)

structure = Structure([x,y], cell, sg)

hamiltonian = Hamiltonian(structure, [HeisenbergExchange(x, y, j=1, cell_offset=(0,0,1))])

hamiltonian.print_summary()

hamiltonian.complete_symmetry()