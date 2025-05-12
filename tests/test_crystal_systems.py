from pyspinw.symmetry.system import crystal_systems

def test_all_systems_present():
    """ Test that all the Bravais lattices are there are correct """

    expected_names = sorted([
        "aP",
        "mP", "mS",
        "oP", "oS", "oI", "oF",
        "tP", "tI",
        "hR",
        "hP",
        "cP", "cI", "cF"
    ])

    assert len(expected_names) == 14 # Sanity check

    names = []
    for crystal_system in crystal_systems:
        names += [crystal_system.letter + bravais.letter for bravais in crystal_system.bravais_options.bravias]

    names.sort()

    for expected, calculated in zip(expected_names, names):
        assert expected == calculated

