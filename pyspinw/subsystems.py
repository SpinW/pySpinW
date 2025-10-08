""" Identification of subsystems of a collection of sites and couplings """
from collections import defaultdict

import numpy as np

from coupling import Coupling
from site import LatticeSite, ImpliedLatticeSite


def find_components(sites: list[LatticeSite], couplings: list[Coupling]):
    """
    Find components (maximal disjoint subsystems) of the given system, or
    in other words, group the sites by being coupled to each other

    Important: Assumes that all the sites given are not parented to sites not in the list"""

    unique_sites = [site for site in sites if not isinstance(site, ImpliedLatticeSite)]
    n_unique_sites = len(unique_sites)
    unique_id_to_index = {unique_sites[index]._unique_id: index for index in range(n_unique_sites)}

    unique_id_to_sites = defaultdict(list[LatticeSite])
    for site in sites:
        unique_id_to_sites[site.parent_site._unique_id].append(site)

    