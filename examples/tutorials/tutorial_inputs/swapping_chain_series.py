""" Antiferromagnetic chain example """

from multiprocessing.spawn import freeze_support

from pyspinw.exchange import HeisenbergExchange
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.interface import generate_exchanges
from pyspinw.path import Path
from pyspinw.site import Site
from pyspinw.symmetry.supercell import SummationSupercell, CommensuratePropagationVector
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structure import Structure
import sys

unit_cell = UnitCell(3, 8, 8)

x = Site(0, 0, 0, 0, 1, 0, name="X")
y = Site(0.5, 0, 0, 0, 1, 0, name="Y")

j1 = HeisenbergExchange(x, y, j=1, cell_offset=(0, 0, 0), name="J1")
j2 = HeisenbergExchange(y, x, j=1, cell_offset=(0, 1, 0), name="J2")

sites = [x, y]
exchanges = [j1, j2]

s = Structure(sites, unit_cell)

hamiltonian = Hamiltonian(s, exchanges)

path = Path([[0,0,0], [0,1,0]])

# One option

# parameterized_hamiltonian = hamiltonian.parameterize(
#     ("J", "j"),
#     find_ground_state_with={"fixed": [x], "verbose": False})


# Another option

# parameterized_hamiltonian = hamiltonian.parameterize(
#     ("J.j"),
#     find_ground_state_with={"fixed": [x], "verbose": False})


# Yet another option

# parameterized_hamiltonian = hamiltonian.parameterize(
#     ((j1, j2), "j"),
#     find_ground_state_with={"fixed": [x], "verbose": False})

# Last option, we'll use this


parameterized_hamiltonian = hamiltonian.parameterize(
    [(j1, "j"), (j2, "j")],
    find_ground_state_with={"fixed": [x], "verbose": False})


print(parameterized_hamiltonian)

j_values = [-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5]

for j in j_values:
    parameterized_hamiltonian(j).print_summary()

parameterized_hamiltonian.energy_plot(j_values, path)
