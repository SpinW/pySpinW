""" Identification of subsystems of a collection of sites and couplings """

from pyspinw.coupling import Coupling
from pyspinw.site import LatticeSite, ImpliedLatticeSite


class Node:
    """ Graph model for depth first search """

    def __init__(self):
        self.connected: list["Node"] = []
        self.group = -1

    def propagate(self, group) -> bool:
        """ Search the graph, only follow if not already visited (self.group>=0)"""
        if self.group == -1:
            self.group = group
            for node in self.connected:
                node.propagate(group)

            return True
        else:
            return False

def find_components(sites: list[LatticeSite], couplings: list[Coupling]) \
        -> list[tuple[list[LatticeSite], list[Coupling]]]:
    """Find components (maximal disjoint subsystems) of the given system, or
    in other words, group the sites by being coupled to each other

    Important: Assumes that all the sites in the couplings are in the list of sites
    """
    unique_sites = [site for site in sites if not isinstance(site, ImpliedLatticeSite)]
    n_unique_sites = len(unique_sites)
    unique_id_to_index = {unique_sites[index]._unique_id: index for index in range(n_unique_sites)}

    # Set up graph
    nodes = [Node() for _ in unique_sites]

    for coupling in couplings:
        a = unique_id_to_index[coupling.site_1.parent_site._unique_id]
        b = unique_id_to_index[coupling.site_2.parent_site._unique_id]

        nodes[a].connected.append(nodes[b])
        nodes[b].connected.append(nodes[a])

    # Depth first graph search
    group = 0
    for node in nodes:
        if node.propagate(group):
            group += 1

    # collect together into groups
    site_groups = [[] for _ in range(group)]
    coupling_groups = [[] for _ in range(group)]

    # Collect sites together
    for site in sites:
        node_index = unique_id_to_index[site.parent_site._unique_id]
        group_index = nodes[node_index].group
        site_groups[group_index].append(site)

    # Collect couplings together
    for coupling in couplings:
        node_index = unique_id_to_index[coupling.site_1.parent_site._unique_id] # Group should be the same for both sites
        group_index = nodes[node_index].group
        coupling_groups[group_index].append(coupling)

    # groups as list of tuples
    return [pair for pair in zip(site_groups, coupling_groups)]
