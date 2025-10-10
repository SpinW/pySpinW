""" Tests for the decomposition into subsystems"""

from pyspinw.site import LatticeSite, ImpliedLatticeSite
from pyspinw.coupling import HeisenbergCoupling

from pyspinw.subsystems import find_components


def test_no_couplings_no_parenting():
    """ Check completely unconnected case returns singletons """
    sites = [LatticeSite(0,0,0),
             LatticeSite(0,0,0),
             LatticeSite(0,0,0)]

    components = find_components(sites, [])

    assert len(components) == len(sites), "Should have same number of components as sites"
    assert all([len(component_sites) == 1 for component_sites, _ in components]), "Components should have one site"

def test_no_couplings_some_parenting():
    """ Check that parented sites get put in same group"""
    # Use the i component to identify two groups
    sites = [LatticeSite(0,0,0),
             LatticeSite(1,0,0)]

    sites += [ImpliedLatticeSite(parent, parent.i, 0, 0) for parent in sites]

    components = find_components(sites, [])

    assert len(components) == 2, "Should have two components"
    assert all([len(component_sites) == 2 for component_sites, _ in components]), "Each component should have two sites"

    # Check the sites are the right sites
    for component_sites, _ in components:
        i = component_sites[0].i
        for site in component_sites:
            assert site.i == i, "Sites within component should be linked by parentage"


def test_two_systems():
    """ Two subsystems defined by couplings with no parenting going on"""
    # expected group encoded in j
    # site within group encoded in i
    sites = [LatticeSite(0,0,0),
             LatticeSite(1,0,0),
             LatticeSite(0,1,0),
             LatticeSite(1,1,0)]

    couplings = [
        HeisenbergCoupling(sites[0], sites[1], 1),
        HeisenbergCoupling(sites[2], sites[3], 1)]

    components = find_components(sites, couplings)

    assert len(components) == 2, "There should be two components"

    for component_sites, component_couplings in components:
        expected_group_index = component_sites[0].j # expected group encoded in j

        assert len(component_sites) == 2, "Each group should contain two sites"
        assert len(component_couplings) == 1, "Each group should contain one coupling"

        for site in component_sites:
            assert site.j == expected_group_index, "Site should be in expected group"

        for coupling in component_couplings:
            assert coupling.site_1.j == expected_group_index, "Coupling should refer to site in expected group"
            assert coupling.site_2.j == expected_group_index, "Coupling should refer to site in expected group"


def test_parenting_connected_two_systems():
    """ Same as the two system test, except the non-coupled sites are connected by parentage,
    therefore there should only be a single group
    """
    sites = [LatticeSite(0, 0, 0),
             LatticeSite(1, 0, 0)]

    sites += [ImpliedLatticeSite(parent, parent.i, 0, 0) for parent in sites]

    couplings = [
        HeisenbergCoupling(sites[0], sites[1], 1),
        HeisenbergCoupling(sites[2], sites[3], 1)]

    components = find_components(sites, couplings)

    assert len(components) == 1, "There should only be one subsystem"

    component_sites, component_couplings = components[0]

    assert len(component_sites) == 4, "The single component should have 4 sites"
    assert len(component_couplings) == 2, "The single component should have 2 couplings"


