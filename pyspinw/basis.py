""" Functions for working with different basis vectors"""

import numpy as np
from pyspinw.dimensionality import dimensionality_check
@dimensionality_check(vectors=(-1, 3))
def find_aligned_basis(vectors: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """ Find a set of basis vectors aligned with the input vectors"""

    return (vectors, vectors, vectors)
