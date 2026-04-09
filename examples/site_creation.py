""" Example site creation """

from pyspinw.interface import generate_sites


s = generate_sites(
    positions=[[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]],
    spins=[[0., 0., i] for i in range(1, 5)],
    names=['X', 'Y', 'Z', 'W'])

x, y, z, w = tuple(s)

print(x)
print(y)
print(z)
print(w)

# Specify supercell spins as a tensor
s = generate_sites(
    positions=[[0, 0, 0], [0.5, 0, 0], [0, 0.5, 0], [0, 0, 0.5]],
    spins=[[[0., 0., i]] for i in range(1, 5)],
    names=['X', 'Y', 'Z', 'W'])

x, y, z, w = tuple(s)

print(x)
print(y)
print(z)
print(w)