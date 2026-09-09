from pyspinw import *

p1 = spacegroup("p1")

p4mmm = spacegroup("p4/mmm")

x = LatticeSite(0.5,0.5,0.5)

p1_ops = p1.operations_between_sites(x, x)
print(p1_ops)

for op in p4mmm.operations_between_sites(x, x):
    print(op)
print()

# Pick a symmetry operation, and make a pair that satisfies it
# 2: -y, x, z

a = LatticeSite(0.2, 0.3, 0.4)
b = LatticeSite(0.7, 0.2, 0.4)

print(p4mmm.operations[2])
print(p4mmm.operations[2]([a.ijk]))

print(p4mmm.operations_between_sites(a, b))

# what about a high symmetry point with the same operation

print()

for op in p4mmm.operations_between_sites(
        LatticeSite(0.5, 0.5, 0.1),
        LatticeSite(0.5,0.5,0.9)):

    print(op)


