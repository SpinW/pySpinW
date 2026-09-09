from pyspinw import *

one_spin_unit_cell = UnitCell(1,1,1)
one_spin_site = LatticeSite(0.5, 0.5, 0.5, sz=1, name="S")

supercell = rotation_supercell(directions=[(0.5, 0, 0)], axes=[(0, 1, 0)])
print(supercell.scaling)

structure = Structure([one_spin_site], one_spin_unit_cell, supercell=supercell)

# We only need to specify one exchange in this case, as it is implicit in the supercell description that
#  exchanges are the same in each repetition. Again, $j=1$ for an antiferromagnet.

one_site_exchange = HeisenbergExchange(one_spin_site, one_spin_site, cell_offset=(1,0,0), j=1)
hamiltonian = Hamiltonian(structure, [one_site_exchange])

hamiltonian.print_summary()
print(hamiltonian.structure.supercell.scaling)

# View it
## skip
view(hamiltonian)
## image: one_site_hamiltonian.png
## snapshot(two_site_structure, "one_site_hamiltonian.png", view_point=(0,-5,-0.5))
