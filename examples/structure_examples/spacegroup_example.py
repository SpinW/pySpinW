from pyspinw import *

structure =Structure(
    [LatticeSite(0, 0, 0, 0, 0, 0)],
    UnitCell(1, 1, 1),
    spacegroup("Im-3m"))

view(structure)