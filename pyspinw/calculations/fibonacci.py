""" Routines for calculating a Fibonacci distribution on a sphere"""
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
    def fibonacci_series(n_numbers):
        """ Generate a Fibonacci series of length `n_numbers`"""
        a, b = 1, 1
        output = []
        for i in range(n_numbers):
            output.append(a)
            a, b = b, a + b

        return output

    @staticmethod
    def fibonacci_to_index(fibonacci_number: float):
        """ Get the first index of the fibonacci series corresponding greater than a given number
        """

        if fibonacci_number < 1:
            raise ValueError("Fibonacci numbers are all bigger than one")

        if fibonacci_number == 1:
            return 1

        # Start with lower bound based on first term of exact formula, as the other term has magnitude less than one
        # the correct index will be either i or i+1.
        i = int(np.log(fibonacci_number * _root_5)/np.log(_golden_ratio) - 1)

        if Fibonacci.fibonacci_by_index(i) < fibonacci_number:
            return i + 1
        else:
            return i

    @staticmethod
    def fibonacci_by_index(i):
        """ Get the ith term of the Fibonacci series"""

        # Two options, Binet's formula or evaluate series,
        # Binets formula gives errors pretty quickly, so might
        # as well use a series

        return Fibonacci.fibonacci_series(i+1)[-1]

    @staticmethod
    def straddling_fibonaccis(n):
        "Get the Fibonacci greater than or equal to n, and the one before "
        last = Fibonacci.fibonacci_to_index(n)

        series = Fibonacci.fibonacci_series(last+1)

        return series[-2], series[-1]

    @staticmethod
    def _phi_z_and_weights(n_points_requested):
        """ Calculate the points and weights for a given requested points"""

        # Notation used in paper, see class docstring

        F_prime, F = Fibonacci.straddling_fibonaccis(n_points_requested)

        delta_z = 2 / F

        j = np.arange(F) # Do we want an extra point?
        z = j*delta_z - 1

        phi = (2 * np.pi * F_prime / F) * j


        weights = 2*(np.pi*delta_z) * (1 + np.cos(np.pi * z))


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
    print(Fibonacci.fibonacci_series(10))
    print([Fibonacci.fibonacci_by_index(i) for i in range(10)])

    for number in Fibonacci.fibonacci_series(10):
        print(Fibonacci.fibonacci_to_index(number))

    for i in range(1, 100):
        print(i, Fibonacci.straddling_fibonaccis(i))

    print(Fibonacci.fibonacci_series(101)[-1])
    print(Fibonacci.fibonacci_by_index(100))


    import matplotlib.pyplot as plt

    phi, z, weights = Fibonacci._phi_z_and_weights(100)
    plt.scatter(phi % (2*np.pi), z)
    plt.show()

