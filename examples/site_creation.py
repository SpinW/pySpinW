""" Example site creation """

from pyspinw.interface import sites


s = sites(
    positions=[[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]],
    moments=[[0., 0., i] for i in range(1, 5)],
    names=['X', 'Y', 'Z', 'W'])

x, y, z, w = tuple(s)

print(x)
print(y)
print(z)
print(w)

# Specify supercell moments as a tensor
s = sites(
    positions=[[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]],
    moments=[[[0., 0., i]] for i in range(1, 5)],
    names=['X', 'Y', 'Z', 'W'])

x, y, z, w = tuple(s)

print(x)
print(y)
print(z)
print(w)