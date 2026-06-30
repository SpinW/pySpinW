"""Generally helpful functions that don't obviously live anywhere else in particular"""
import functools
from collections import defaultdict
from typing import TypeVar

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
                    # Opposite spins

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


def rotation_from_z(target_vector):
    """ Rotation matrix from (0,0,1) to the target vector direction """
    mag_sq = np.sum(target_vector**2)

    if mag_sq < 1e-9:
        return np.eye(3)

    v = target_vector / np.sqrt(mag_sq)

    x,y,z = v

    if z < 1e-9 - 1:
        return np.array([[-1,0,0],[0,1,0],[0,0,-1]])

    return np.array([
        [1-x**2 / (1+z), -x*y / (1+z), x],
        [-x*y / (1+z), 1 - y**2 / (1+z), y],
        [-x, -y, z]
    ])

def is_diagonal(m: np.ndarray):
    """ Check whether an array is diagonal """
    return np.all(m == np.diag(np.diagonal(m)))


T = TypeVar("T")
class IncrementalApproximateHistogram[T]:
    """ Histogramming method for finding nearest neighbour distances in lattices """

    def __init__(self, tolerance=1e-9):
        self._tolerance = tolerance
        self._groups = defaultdict(list)

    def add(self, distance, entry: T):
        """ Add an entry to the collection"""
        for key in self._groups:
            if abs(key - distance) < self._tolerance:
                self._groups[key].append(entry)

        self._groups[distance].append(entry)

    def groups(self):
        """ Groups in order of distance"""
        ordered_keys = sorted(self._groups.keys())
        return [self.groups[key] for key in ordered_keys]
    def group_sizes(self):
        """ Size of each group"""
        return [len(group) for group in self.groups()]
    def n_groups(self):
        """ Number of groups at different distances """
        return len(self._groups)


@functools.lru_cache(maxsize=10)
def cell_shell(order: int) -> list[tuple[int, int, int]]:
    """ Create cell offsets for shells of unit cells

    i.e.


      (0)      (1)      (2)      (3)
    -------  -------  -------  #######
    -------  -------  -#####-  #-----#
    -------  --###--  -#---#-  #-----#
    ---#---  --#-#--  -#---#-  #-----#
    -------  --###--  -#---#-  #-----#
    -------  -------  -#####-  #-----#
    -------  -------  -------  #######

    """
    # This is basically all the integer vectors with a L_infty distance of a given integer
    # L_infty = max(x,y,z)

    out = []
    for x in range(-order, order+1):
        for y in range(-order, order+1):
            for z in range(-order, order+1):
                if max([abs(x), abs(y), abs(z)]) == order:
                    out.append((x, y, z))

    return out
