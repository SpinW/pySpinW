""" Replicates tutorial 20 """

import numpy as np
from pyspinw import *

symmetries = '-z, y+3/4, x+3/4; z+3/4, -y, x+3/4; z+3/4, y+3/4, -x; y+3/4, x+3/4, -z; x+3/4, -z, y+3/4; -z, x+3/4, y+3/4'

sg = get_spacegroup(symmetries) # Gives spacegroup F 4_1/d -3 2/m, second setting, hall=526
print(sg)

a = 10.0307
unit_cell = UnitCell(a, a, a)

g_tensor = 4.32*np.eye(3) - 0.84

print(g_tensor)

yb = Site(1 / 2, 1 / 2, 1 / 2, 0, 0, 1 / 2, name="Yb 3+", g=g_tensor)
ti = Site(0, 0, 0, name="Ti 4+")
o1 = Site(0.3318, 1 / 8, 1 / 8, name="O 2-")
o2 = Site(3 / 8, 3 / 8, 3 / 8, name="O 2- ")

sites = [yb, ti, o1, o2]

structure = Structure(sites, unit_cell, sg)

structure.print_summary()

yb_neighbour, offset = structure.neighbours(yb, element="Yb")[-1]

hamiltonian = Hamiltonian(structure, [], [])

paths = map(Path, [
    [[-0.5, -0.5, -0.5], [2, 2, 2]],
    [[1, 1, -2], [1, 1, 1.5]],
    [[2, 2, -2], [2, 2, 1.5]],
    [[-0.5, -0.5, 0], [2.5, 2.5, 0]],
    [[0, 0, 1], [2.3, 2.3, 1]]
])

field_direction = np.array([1, -1, 0], dtype=float)
field_direction /= np.sqrt(np.sum(field_direction**2))
for field_magnitude in [2, 5]:
    field = field_direction * field_magnitude

    optimised = hamiltonian.ground_state(field=field)

    view(optimised)