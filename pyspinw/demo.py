""" Demonstrations / packaging tests """

def run_demos():
    """ Run through general tasks that should test things are installed correctly"""
    demo_viewer()
    demo_chains()

def demo_chains():
    """ Antiferromagnetic chain example """

    from pyspinw import (UnitCell, LatticeSite, HeisenbergCoupling, Structure,
                         Hamiltonian, Path)

    unit_cell = UnitCell(3, 8, 8)

    x = LatticeSite(0, 0, 0, 0, 1, 0, name="X")
    y = LatticeSite(0.5, 0, 0, 0, 1, 0, name="Y")

    j1 = HeisenbergCoupling(x, y, j=1, cell_offset=(0, 0, 0), name="J1")
    j2 = HeisenbergCoupling(y, x, j=1, cell_offset=(0, 1, 0), name="J2")

    sites = [x, y]
    exchanges = [j1, j2]

    s = Structure(sites, unit_cell)

    hamiltonian = Hamiltonian(s, exchanges)

    path = Path([[0, 0, 0], [0, 1, 0]])

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

def demo_viewer():
    """ Show the viewer with an example """
    from pyspinw import (
        UnitCell, LatticeSite, Structure, TrivialSupercell,
        generate_exchanges, HeisenbergCoupling, Hamiltonian, view)

    unit_cell = UnitCell(1,1,1, gamma=60)

    x = LatticeSite(0, 0, 0, 0, 0, 1, name="X")
    y = LatticeSite(0.5, 0, 0, 0, 0, 1, name="Y")
    z = LatticeSite(0, 0.5, 0, 0, 0, 1, name="Z")
    w = LatticeSite(0.5, 0.5, 0, 0, 0, 0, name="W")

    sites = [x, y, z, w]

    s = Structure(sites, unit_cell=unit_cell, supercell=TrivialSupercell(scaling=(3,3,1)))

    exchanges = generate_exchanges(sites=[x, y, z],
                                   unit_cell=unit_cell,
                                   max_distance=0.6,
                                   coupling_type=HeisenbergCoupling,
                                   j=-1)

    hamiltonian = Hamiltonian(s, exchanges)

    hamiltonian.print_summary()

    view(hamiltonian)

