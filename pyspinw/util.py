"""Generally helpful functions that don't obviously live anywhere else in particular"""

import numpy as np

from pyspinw.checks import check_sizes
from pyspinw.site import LatticeSite
from pyspinw.tolerances import tolerances


@check_sizes(v=(3,), force_numpy=True)
def triple_product_matrix(v: np.ndarray):
    """Find the matrix, M, such that for all vectors X, Y

    X^T M Y = V . (X x Y)
    """
    x, y, z = v

    return np.array([
        [ 0,  z, -y],
        [-z,  0,  x],
        [ y, -x,  0]])


def demo_triple_product_matrix():
    """Show an example of making a matrix that does the triple product"""
    v = [1, 2, 3]
    m = triple_product_matrix(v)
    print(m)

def problematic_sites(sites: list[LatticeSite], implied_sites: list[LatticeSite], site_to_implied_site: list[list[int]]):
    """ Find sites that are contradictory in terms of magnetism (i.e. they have an image with same position but
    negated spin)

    :param sites: Asymmetric cell sites
    :param implied_sites: Sites inferred by symmetry
    :param site_to_implied_site: mapping from site index to implied site index

    :returns: Indices of sites with this problem

    """
    bad_sites = []

    for i, site in enumerate(sites):
        for j in site_to_implied_site[i]:
            implied_site = implied_sites[j]
            if np.all(np.abs(site.ijk - implied_site.ijk) < tolerances.SAME_SITE_ABS_TOL):
                # Same position

                if np.all(np.abs(site.m + implied_site.m) < tolerances.SAME_SITE_ABS_TOL):
                    # Opposite moments

                    bad_sites.append(i)
                    break

    return bad_sites


if __name__ == "__main__":
    demo_triple_product_matrix()
