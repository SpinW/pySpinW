import numpy as np

from pyspinw.checks import check_sizes
from pyspinw.supercell import CommensurateSupercell, CommensuratePropagationVector


def _rotation_transformation_supercell(propagation_vectors: list[CommensuratePropagationVector],
                                       axis: np.ndarray,
                                       orthogonal: bool):

    rotations = 

def _rotation_summation_supercell(k: np.ndarray, axis: np.ndarray, initial_moment: bool, orthogonal: bool):
    pass

@check_sizes(k=(3,))
def rotation_supercell(k: np.ndarray,
                       axis: np.ndarray | None =None,
                       initial_moment: np.ndarray | None = None,
                       orthogonal: bool=False):

    if axis is None:
        axis = np.array([0,0,1])

    if orthogonal:
        k_vectors = [CommensuratePropagationVector(v) for v in np.diag(k)]
    else:
        k_vectors = [CommensuratePropagationVector(k)]

    if initial_moment is None:
        return _rotation_transformation_supercell(k_vectors, axis, orthogonal)
    else:
        return _rotation_summation_supercell(k_vectors, axis, initial_moment)


if __name__ == "__main__":
    pass