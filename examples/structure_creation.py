from pyspinw.interface import spacegroup
from pyspinw.site import LatticeSite
from pyspinw.symmetry.supercell import TrivialSupercell
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structures import Structure

from pyspinw.debug_plot import debug_plot

print("Cubic high-symmetry group")

unit_cell = UnitCell(1,1,1)
group = spacegroup("i4/m -3 2/m")

s = Structure([
    LatticeSite(0.1, 0.1, 0.1),
    LatticeSite(0.2, 0.0, 0.2)],
    unit_cell=unit_cell,
    spacegroup=group,
    supercell=TrivialSupercell())

for site in s.sites:
    print(site)

debug_plot(s, [], show=False)

print("P1 with supercell")

unit_cell = UnitCell(1,1,1)
group = spacegroup("p1")

s = Structure([
    LatticeSite(0.1, 0.1, 0.1),
    LatticeSite(0.2, 0.0, 0.2)],
    unit_cell=unit_cell,
    spacegroup=group,
    supercell=TrivialSupercell(scaling=(3,3,3)))

for site in s.full_structure_site_list():
    print(site)


debug_plot(s, [])
