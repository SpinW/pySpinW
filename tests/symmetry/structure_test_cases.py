from pyspinw import *

def orthogonal_unit_cell_transform_p4mm():
    sg = spacegroup("p4mm")
    sites = [LatticeSite(0.2, 0.2, 0.5, 0,0,1, name="site")]
    unit_cell = UnitCell(1,1,2)

    return Structure(sites, unit_cell=unit_cell, spacegroup=sg)

def non_orthogonal_unit_cell_transform_p6mm():
    sg = spacegroup("p6mm")
    sites = [LatticeSite(0.2, 0.2, 0.5, 0,0,1, name="site")]
    unit_cell = UnitCell(1,1,2, gamma=120)

    return Structure(sites, unit_cell=unit_cell, spacegroup=sg)




cases = [orthogonal_unit_cell_transform_p4mm, non_orthogonal_unit_cell_transform_p6mm]

if __name__ == "__main__":
    for case in cases:
        view(case())