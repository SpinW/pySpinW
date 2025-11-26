"""Generally helpful functions that don't obviously live anywhere else in particular"""

import numpy as np
from numpy._typing import ArrayLike

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

def problematic_sites(sites: list[LatticeSite],
                      implied_sites: list[LatticeSite],
                      site_to_implied_site: list[list[int]]):
    """ Find sites that are contradictory in terms of magnetism

    These sites have an image with same position but negated spin

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

@check_sizes(axis=(3,), force_numpy=True)
def rotation_matrix(angle, axis):
    """ Create a rotation matrix """
    mag = np.sqrt(np.sum(axis**2))

    if mag == 0:
        raise ValueError("Rotation matrix cannot be made for axis (0,0,0) ")

    axis = axis/mag  # Don't use /= because of dtype error possibility

    c = np.cos(angle)
    s = np.sin(angle)

    part1 = (axis.reshape(-1, 1) * axis.reshape(1, -1)) * (1-c)
    part2 = c * np.eye(3)
    part3 = s * triple_product_matrix(-axis)

    return part1 + part2 + part3

def connected_components(adjacency_matrix: np.ndarray) -> list[list[int]]:
    """ Get the connected components of a graph specied by an adjacency matrix

    :param adjacency_matrix: n-by-n numpy array of dtype bool representing the adjacency matrix
    :returns: list of subgraphs, themselves lists of indices for the adjacency matrix
    """

    assert len(adjacency_matrix.shape) == 2
    assert adjacency_matrix.shape[0] == adjacency_matrix.shape[1]

    n = adjacency_matrix.shape[0]
    visited = np.zeros((n, ), dtype=bool)
    components = []

    def search(current_node, component):
        visited[current_node] = True
        component.append(current_node)
        for i in range(n):
            if adjacency_matrix[current_node][i] and not visited[i]:
                search(i, component)

    for i in range(n):
        if not visited[i]:
            component = []
            search(i, component)
            components.append(component)

    return components

def arraylike_equality(array_1: ArrayLike, array_2: ArrayLike, abs_tol=1e-8):
    """ Check for approximate equality of arrays that may or may not have the same shape """
    array_1 = np.array(array_1)
    array_2 = np.array(array_2)

    if array_1.shape != array_2.shape:
        return False

    return np.all(np.abs(array_1 - array_2) < abs_tol)

if __name__ == "__main__":
    demo_triple_product_matrix()
