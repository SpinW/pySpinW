""" Routines for calculating a Fibonacci distribution on a sphere

See: J H Hannay and J F Nye 2004 J. Phys. A: Math. Gen. 37 11591
DOI: 10.1088/0305-4470/37/48/005

"""


import numpy as np

_root_5 = np.sqrt(5)
_golden_ratio = 0.5*(_root_5 + 1)
_golden_ratio_inv_neg = -0.5*(_root_5 - 1)


class Fibonacci:
    """ Routines for calculating a Fibonacci distribution on a sphere

    See: J H Hannay and J F Nye 2004 J. Phys. A: Math. Gen. 37 11591
    DOI: 10.1088/0305-4470/37/48/005

    """

    @staticmethod
    def straddling_fibonaccis(n):
        "Get the Fibonacci greater than or equal to n, and the one before "

        a = 1
        b = 1

        while b < n:
            a, b = b, a + b

        return a, b

    @staticmethod
    def _phi_z_and_weights(n_points_requested):
        """ Calculate the points and weights for a given requested points"""
        # Notation used in paper, see class docstring

        F_prime, F = Fibonacci.straddling_fibonaccis(n_points_requested)

        delta_z = 2 / F

        j = np.arange(F) # Do we want an extra point?
        z = j*delta_z - 1

        phi = (2 * np.pi * F_prime / F) * j


        # weights = 2*(np.pi*delta_z) * (1 + np.cos(np.pi * z))
        weights = np.ones((F, )) * (4*np.pi / F)

        return phi, z, weights

    @staticmethod
    def points_and_weights(n_points_requested):
        """ Get points on sphere and weights """
        # Get points and weights in phi, z coordinates
        phi, z, weights = Fibonacci._phi_z_and_weights(n_points_requested)

        # Convert phi and z into 3D points

        points = np.zeros((phi.shape[0], 3))

        r = np.sqrt(1-z**2)

        points[:, 0] = r*np.sin(phi)
        points[:, 1] = r*np.cos(phi)
        points[:, 2] = z

        # Weights *should* be the same

        return points, weights


if __name__ == "__main__":

    for i in range(1, 100):
        print(i, Fibonacci.straddling_fibonaccis(i))


    import matplotlib.pyplot as plt

    phi, z, weights = Fibonacci._phi_z_and_weights(100)
    plt.scatter(phi % (2*np.pi), z)
    plt.show()

