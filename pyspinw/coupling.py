from idlelib.squeezer import Squeezer

import numpy as np
from numpy.random.mtrand import Sequence

from pyspinw._base import MagneticStructure


class Coupling:
    """ Coupling between pairs of atoms """

    def __init__(self,
                 structure: MagneticStructure,
                 site_1_identifier: str, # TODO, how will we identify atoms
                 site_2_identifier: str,
                 coupling_matrix: int | float | Sequence[float] | Sequence[Sequence[float]] | np.ndarray):

        # If it's a number, convert to a numpy array
        if isinstance(coupling_matrix, (int, float, Sequence)):
            coupling_matrix = np.array([coupling_matrix])

        if isinstance(coupling_matrix, np.ndarray):
            match coupling_matrix.shape:
                case (1,) | ():
                    # A single number, same coupling in each direction, i.e. J is a scalar constant
                    coupling_matrix = np.eye(3, dtype=float) * coupling_matrix

                case (3,):
                    # Diagonal elements of J specified
                    coupling_matrix = np.diag(v=coupling_matrix)

                case (3, 3):
                    # Full matrix, do nothing
                    pass

                case _:
                    raise ValueError("Expected coupling_matrix to be a scalar, 3-vector, or 3-by-3 matrix")

        else:
            raise TypeError("Expected coupling_matrix to be an int, float or numpy array")