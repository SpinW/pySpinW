""" Antiferromagnetic chain example """


from pyspinw.coupling import HeisenbergCoupling
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structures import Structure

unit_cell = UnitCell(3, 8, 8)

x = LatticeSite(0, 0, 0, 0, 1, 0, name="X")
y = LatticeSite(0.5, 0, 0, 0, 1, 0, name="Y")

j1 = HeisenbergCoupling(x, y, j=1, cell_offset=(0,0,0), name="J1")
j2 = HeisenbergCoupling(y, x, j=1, cell_offset=(0,1,0), name="J2")

sites = [x, y]
exchanges = [j1, j2]

s = Structure(sites, unit_cell)

hamiltonian = Hamiltonian(s, exchanges)

parameterized_hamiltonian = hamiltonian.parameterize(
    ("J", "j"),
    find_ground_state_with={"fixed": [x], "verbose": False})

