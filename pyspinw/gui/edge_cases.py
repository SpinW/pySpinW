import numpy as np

def _hypercube(n):
    """ Inner loop for `hypercube`"""
    if n <= 0:
        return [[]]

    else:
        smaller = _hypercube(n-1)
        return [x + [0] for x in smaller] + [x + [1] for x in smaller]

def hypercube(n):
    """ Generate points in an n-dimensional unit hypercube"""
    return [np.array(x) for x in _hypercube(n)]

# rules for replacement
replacements = {i: hypercube(i) for i in range(4)}


def add_extra_edge_lines(point_pairs: list[tuple[np.ndarray, np.ndarray]], tol=1e-12):
    """ When a line lives on the x/y/z=0 face/edges of the unit cube, make extra ones

    * If a line lives on a face an not an edge, i.e. one component in both points are zero, we replace
       two copies, one where the component is zero, the other where it is one.
    * If a line lives on an edge, i.e. 2 two of the components in both points are zero, we need to make 4 copies
       with the two components being [0,0], [0,1], [1,0] and [0,0]
    * Three being zero should not happen, because the line will have length zero
    """

    new_pair_list = []
    for a, b in point_pairs:
        a_zeros = np.isclose(a, 0, atol=tol)
        b_zeros = np.isclose(b, 0, atol=tol)

        both_zeros = np.logical_and(a_zeros, b_zeros)

        n_zeros = np.sum(both_zeros)

        for replacement in replacements[n_zeros]:
            a_copy = a.copy()
            b_copy = b.copy()
            a_copy[both_zeros] = replacement
            b_copy[both_zeros] = replacement
            new_pair_list.append((a_copy, b_copy))

    return new_pair_list



def add_extra_edge_points(point, tol=1e-12):
    """ If a point has any coordinate of zero, add extra ones at 1 """

    zeros = np.isclose(point, 0, atol=tol)

    n_zeros = int(np.sum(zeros))

    new_points = []
    for replacement in replacements[n_zeros]:
        point_copy = point.copy()
        point_copy[zeros] = replacement

        new_points.append(point_copy)

    return new_points