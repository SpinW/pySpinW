from pyspinw import *
from pyspinw import HeisenbergExchange

# symmetries='-z, y+3/4, x+3/4; z+3/4, -y, x+3/4; z+3/4, y+3/4, -x; y+3/4, x+3/4, -z; x+3/4, -z, y+3/4; -z, x+3/4, y+3/4'
#
# symmetries = '-x, y, z; x, -y, z; y, x, z'
# sg = spacegroup(symmetries)

def heisenberg_case_1():
    """ Test case, cubic unit cell """

    sg = spacegroup("p4mm")

    s = LatticeSite(0.25,0.25,0.25,0,0,1, name="X")

    structure = Structure([s], unit_cell=UnitCell(1,1,1), spacegroup=sg)

    s1 = structure.site_by_name("X")
    s2 = structure.site_by_name("X [3]")

    ex1 = HeisenbergExchange(s1, s2, cell_offset=(0,0,0), j=-1)
    ex2 = HeisenbergExchange(s1, s2, cell_offset=(0,-1,0), j=-1)

    hamiltonian = Hamiltonian(structure, [ex1, ex2])

    return hamiltonian

def heisenberg_case_2():
    """ Test case, rhombic unit cell """

    sg = spacegroup("p6mm")

    a = LatticeSite(0, 0, 0, 0, 0, 1, name="a")
    b = LatticeSite(0.5, 0, 0, 0, 1, 0, name="b")
    c = LatticeSite(0, 0.5, 0, 1, 0, 0, name="c")

    structure = Structure([a, b, c], unit_cell=UnitCell(1, 1, 1, gamma=120), spacegroup=sg)


    exchanges = [
        HeisenbergExchange(a, b, cell_offset=(0, 0, 0), j=-1, name="ab"),
        HeisenbergExchange(b, c, cell_offset=(0, 0, 0), j=-1, name="bc"),
        HeisenbergExchange(c, a, cell_offset=(0, 0, 0), j=-1, name="ca"),
    ]

    hamiltonian = Hamiltonian(structure, exchanges)

    return hamiltonian


heisenberg_cases = [
    heisenberg_case_1(),
    heisenberg_case_2(),
]

if __name__ == "__main__":
    for case in [heisenberg_case_2]:
    # for case in [case_1, case_2]:

        h = case()

        view(h)
        view(h.symmetry_filled())