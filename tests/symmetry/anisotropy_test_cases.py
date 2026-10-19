""" Test cases for anisotropy symmetry stuff """
import numpy as np

from pyspinw import UnitCell, spacegroup, LatticeSite, Structure, Hamiltonian, Anisotropy
from pyspinw.interface import view
from pyspinw.symmetry.symbolic import evaluate

def anisotropy_test_system():
    """ Test system with only symmetry unrelated sites, covering all Wyckoff positions """

    unit_cell = UnitCell(1, 1, 2)
    sg = spacegroup("P4mm")

    anisotropy_matrix_values = {letter: number for number, letter in enumerate("abcdef", start=1)}

    # Pick sites at wyckoff positions in particular, none should be symmetry related
    sites = [
        LatticeSite(0, 0, 0.5, name="a/0.5"),  # a / 0.5
        LatticeSite(0.5, 0.5, 0.5, name="b/0.5"),  # b / 0.5
        LatticeSite(0.5, 0.5, 0.3, name="b/0.3"),  # b / 0.3
        LatticeSite(0.5, 0, 0.5, name="c/0.5"),  # c / 0.5
        LatticeSite(0.5, 0, 0.3, name="c/0.3"),  # c / 0.3
        LatticeSite(0.25, 0.25, 0.4, name="d/0.25,0.4"),  # d / 0.25,0.4
        LatticeSite(0.35, 0.35, 0.4, name="d/0.35,0.4"),  # d / 0.35, 0.4
        LatticeSite(0.25, 0, 0.25, name="e/0.25,0.5"),  # e / 0.25,0.25
        LatticeSite(0.15, 0, 0.35, name="e/0.15,0.35"),  # e / 0.15,0.35
        LatticeSite(0.25, 0, 0.35, name="e/0.25,0.35"),  # e / 0.25,0.35
        LatticeSite(0.4, 0.5, 0.1, name="f/0.4,0.1"),  # f / 0.4, 0.1
        LatticeSite(0.7, 0.5, 0.1, name="f/0.7,0.1"),  # f / 0.7, 0.1
        LatticeSite(0.8, 0.5, 0.1, name="f/0.8,0.1"),  # f / 0.8, 0.1
        LatticeSite(0.8, 0.5, 0.2, name="f/0.8,0.2"),  # f / 0.8, 0.2
        LatticeSite(0.1, 0.2, 0.3, name="g1"),  # g
        LatticeSite(0.1, 0.2, 0.35, name="g2"),  # g
        LatticeSite(0.15, 0.2, 0.35, name="g3"),  # g
        LatticeSite(0.4, 0.2, 0.3, name="g4"),  # g
    ]

    for site in sites:
        site.spin = 0, 0, 1

    structure = Structure(sites, unit_cell=unit_cell, spacegroup=sg)

    anisotropies = []
    for site in sites:
        allowable = structure.anisotropy_constraints(site, do_print=False)
        expressions = [s.strip() for s in allowable._matrix_form_strings()[0]]

        v = [evaluate(expression, anisotropy_matrix_values) for expression in expressions]
        matrix = np.array([[v[0], v[1], v[2]],
                           [v[1], v[3], v[4]],
                           [v[2], v[4], v[5]]
                           ], dtype=float)

        anisotropies.append(Anisotropy(site, matrix))

    for anisotropy in anisotropies:
        assert anisotropy.obeys_symmetry(structure), "Anisotropy should obey symmetry"

    return Hamiltonian(structure, [], anisotropies=anisotropies)


if __name__ == "__main__":
    hamiltonian = anisotropy_test_system()
    for anisotropy in hamiltonian.anisotropies:
        print(anisotropy)

    view(hamiltonian.symmetry_filled())