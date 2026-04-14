from pyspinw import *

structure = Structure(
    [LatticeSite(0.5, 0.5, 0.5, 0, 0, 1)],
    UnitCell(1, 1, 1),
    spacegroup("p1"),
    supercell=TiledSupercell(scaling=(3,6,1)))

view(structure)