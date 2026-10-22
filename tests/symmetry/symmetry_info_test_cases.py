from pyspinw import get_spacegroup, UnitCell, Site, Structure
from pyspinw.interface import view


def case_1() -> Structure:
    unit_cell = UnitCell(1,1,2)
    sg = get_spacegroup("P4mm")

    # Pick sites at wyckoff positions in particular
    sites = [
        Site(0, 0, 0.5, name="a/0.5"), # a / 0.5
        Site(0.5, 0.5, 0.5, name="b/0.5"), # b / 0.5
        Site(0.5, 0.5, 0.3, name="b/0.3"), # b / 0.3
        Site(0.5, 0, 0.5, name="c/0.5"), # c / 0.5
        Site(0.5, 0, 0.3, name="c/0.3"), # c / 0.3
        Site(0.25, 0.25, 0.4, name="d/0.25,0.4"), # d / 0.25,0.4
        Site(0.35, 0.35, 0.4, name="d/0.35,0.4"), # d / 0.35, 0.4
        Site(0.25, 0, 0.25, name="e/0.25,0.5"), # e / 0.25,0.25
        Site(0.15, 0, 0.35, name="e/0.15,0.35"), # e / 0.15,0.35
        Site(0.25, 0, 0.35, name="e/0.25,0.35"), # e / 0.25,0.35
        Site(0.4, 0.5, 0.1, name="f/0.4,0.1"), # f / 0.4, 0.1
        Site(0.7, 0.5, 0.1, name="f/0.7,0.1"), # f / 0.7, 0.1
        Site(0.8, 0.5, 0.1, name="f/0.8,0.1"), # f / 0.8, 0.1
        Site(0.8, 0.5, 0.2, name="f/0.8,0.2"), # f / 0.8, 0.2
        Site(0.1, 0.2, 0.3, name="g1"), # g
        Site(0.1, 0.2, 0.35, name="g2"), # g
        Site(0.15, 0.2, 0.35, name="g3"), # g
        Site(0.4, 0.2, 0.3, name="g4"), # g
        ]

    return Structure(sites, unit_cell=unit_cell, spacegroup=sg)

cases = [case_1()]


if __name__ == "__main__":
    s = case_1()
    # view(s)
    #
    # for site_1, site_2 in s.site_pairs():
    #     s.exchange_constraints(site_1, site_2)

    for site in s.sites:
        s.anisotropy_constraints(site)