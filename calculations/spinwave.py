import numpy as np

from pyspinw.checks import check_sizes


@check_sizes(rotations=('n',3,3),
             magnitudes=('n',3),
             q_vectors=('m', 3),
             coupling_r=('p', 3),
             coupling_matrices=('p', 3, 3))
def spinwave_calculation(rotations: np.ndarray,
                         magnitudes: np.ndarray,
                         q_vectors: np.ndarray,
                         coupling_r: np.ndarray,
                         coupling_matrices: np.ndarray):
    """ Main calculation step

    Unlike the main interface it takes indexed arrays, the meaning of the arrays is set elsewhere

    """

    # matrix used for calculating the phase factors, q.r in exp(i q.r) - could be too big, but we'll see
    qr = np.dot(coupling_r, q_vectors)

    phase_factors = np.exp(1j * qr)

    z = rotations[:,:,0] + 1j*rotations[:,:,1] # n-by-3, complex
    eta = rotations[:,:,2] # n-by-3 real

    # Create the A, B, and C matrices

