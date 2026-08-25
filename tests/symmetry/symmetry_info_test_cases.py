from pyspinw import spacegroup, UnitCell, LatticeSite, Structure
from pyspinw.interface import view


def case_1() -> Structure:
    unit_cell = UnitCell(1,1,2)
    sg = spacegroup("P4mm")

    # Pick sites at wyckoff positions in particular
    sites = [
        LatticeSite(0,0,0.5, name="a/0.5"), # a / 0.5
        LatticeSite(0.5, 0.5, 0.5, name="b/0.5"), # b / 0.5
        LatticeSite(0.5, 0.5, 0.3, name="b/0.3"), # b / 0.3
        LatticeSite(0.5,0,0.5, name="c/0.5"), # c / 0.5
        LatticeSite(0.5,0,0.3, name="c/0.3"), # c / 0.3
        LatticeSite(0.25,0.25, 0.4, name="d/0.25,0.4"), # d / 0.25,0.4
        LatticeSite(0.35,0.35, 0.4, name="d/0.35,0.4"), # d / 0.35, 0.4
        LatticeSite(0.25,0, 0.25, name="e/0.25,0.5"), # e / 0.25,0.25
        LatticeSite(0.15,0, 0.35, name="e/0.15,0.35"), # e / 0.15,0.35
        LatticeSite(0.25,0, 0.35, name="e/0.25,0.35"), # e / 0.25,0.35
        LatticeSite(0.4, 0.5, 0.1, name="f/0.4,0.1"), # f / 0.4, 0.1
        LatticeSite(0.7, 0.5, 0.1, name="f/0.7,0.1"), # f / 0.7, 0.1
        LatticeSite(0.8, 0.5, 0.1, name="f/0.8,0.1"), # f / 0.8, 0.1
        LatticeSite(0.8, 0.5, 0.2, name="f/0.8,0.2"), # f / 0.8, 0.2
        LatticeSite(0.1,0.2,0.3, name="g1"), # g
        LatticeSite(0.1,0.2,0.35, name="g2"), # g
        LatticeSite(0.15,0.2,0.35, name="g3"), # g
        LatticeSite(0.4,0.2,0.3, name="g4"), # g
        ]

    return Structure(sites, unit_cell=unit_cell, spacegroup=sg)


if __name__ == "__main__":
    s = case_1()
    # view(s)

    for site_1, site_2 in s.site_pairs():
        s.exchange_constraints(site_1, site_2)