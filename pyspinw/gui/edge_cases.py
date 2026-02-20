import numpy as np


def add_extra_edge_lines(point_pairs: list[tuple[np.ndarray, np.ndarray]], tol=1e-12):
    """ When a line lives on the x/y/z=0 face/edges of the unit cube, make extra ones

    * If a line lives on a face an not an edge, i.e. one component in both points are zero, we replace
       two copies, one where the component is zero, the other where it is one.
    * If a line lives on an edge, i.e. 2 two of the components in both points are zero, we need to make 4 copies
       with the two components being [0,0], [0,1], [1,0] and [0,0]
    * Three being zero should not happen, because the line will have length zero
    """

    # Extras for edges
    subs = [np.array([0.0, 1.0]),
            np.array([1.0, 0.0]),
            np.array([1.0, 1.0])]

    new_pair_list = []
    for a, b in point_pairs:
        a_zeros = np.isclose(a, 0, atol=tol)
        b_zeros = np.isclose(b, 0, atol=tol)

        both_zeros = np.logical_and(a_zeros, b_zeros)

        n_zeros = np.sum(both_zeros)

        new_pair_list.append((a, b)) # We'll always need this
        match n_zeros:
            case 0:
                pass
            case 1:
                a_copy = a.copy()
                b_copy = b.copy()
                a_copy[both_zeros] = 1.0
                b_copy[both_zeros] = 1.0
                new_pair_list.append((a_copy, b_copy))
            case 2:
                for sub in subs:
                    a_copy = a.copy()
                    b_copy = b.copy()
                    a_copy[both_zeros] = sub
                    b_copy[both_zeros] = sub
                    new_pair_list.append((a_copy, b_copy))

            case _:
                raise ValueError("More than 2 shared zeros in pair, either they're the same point"
                                 "or something much worse has happened")

    return new_pair_list