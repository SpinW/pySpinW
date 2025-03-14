import numpy as np

from pyspinw.checks import check_sizes


@check_sizes(rotations=('n',3,3), q_vectors=('m', 3))
def spinwave_calculation(rotations: np.ndarray, q_vectors: np.ndarray):
    """ Main calculation step

    Unlike the main interface it takes indexed arrays, the meaning of the arrays is set elsewhere

    """

    z = rotations[:,:,0] + 1j*rotations[:,:,1] # n-by-3, complex
    eta = rotations[:,:,2] # n-by-3 real

    # Create the A, B, and C matrices
