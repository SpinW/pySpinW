from pyspinw.couplinggroup import CouplingGroup
from pyspinw.interface import spacegroup
from pyspinw.site import LatticeSite
from pyspinw.symmetry.supercell import TrivialSupercell
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structures import Structure

from pyspinw.debug_plot import debug_plot


unit_cell = UnitCell(1,1,1)
group = spacegroup("p1")

s = Structure([LatticeSite(0, 0, 0)],
    unit_cell=unit_cell,
    spacegroup=group,
    supercell=TrivialSupercell(scaling=(3,3,3)))
    