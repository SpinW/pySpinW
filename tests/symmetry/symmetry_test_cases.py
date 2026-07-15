from pyspinw import *
from pyspinw import HeisenbergExchange

# symmetries='-z, y+3/4, x+3/4; z+3/4, -y, x+3/4; z+3/4, y+3/4, -x; y+3/4, x+3/4, -z; x+3/4, -z, y+3/4; -z, x+3/4, y+3/4'
#
# symmetries = '-x, y, z; x, -y, z; y, x, z'
# sg = spacegroup(symmetries)

def case_1():

    sg = spacegroup("p4mm")

    s = LatticeSite(0.25,0.25,0.25,0,0,1, name="X")

    structure = Structure([s], unit_cell=UnitCell(1,1,1), spacegroup=sg)

    s1 = structure.site_by_name("X")
    s2 = structure.site_by_name("X [3]")

    ex1 = HeisenbergExchange(s1, s2, cell_offset=(0,0,0), j=-1)
    ex2 = HeisenbergExchange(s1, s2, cell_offset=(0,-1,0), j=-1)

    hamiltonian = Hamiltonian(structure, [ex1, ex2])

    return hamiltonian

if __name__ == "__main__":
    h = case_1()

    view(h)
    view(h.symmetry_filled())