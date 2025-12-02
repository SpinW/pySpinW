""" Spacegroup and unit cell creation examples """

from pyspinw.interface import spacegroup

group = spacegroup("r-3h")
print(group)
print(group.lattice_system)
print(group.create_unit_cell(a=1, c=2))
print()


group = spacegroup("r-3r")
print(group)
print(group.lattice_system)
print(group.create_unit_cell(1, 99))
print()

group = spacegroup("p3")
print(group)
print(group.lattice_system)
print(group.create_unit_cell(1, 2))
print()

group = spacegroup("x,y,z")
print(group)
print(group.lattice_system)
print(group.create_unit_cell(1,2,3,50,60,80))
print()

# Lots of groups have this, errors out
#group = spacegroup("-x,-y,-z")

group = spacegroup("x,y,z; -x,-y,-z")
print(group)
print(group.lattice_system)
print(group.create_unit_cell(1,2,3,50,60,80))
print()

