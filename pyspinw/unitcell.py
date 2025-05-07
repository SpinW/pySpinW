"""A class to define a unit cell."""
import numpy as np

from pyspinw.checks import check_sizes


class UnitCell:
    """Unit Cell Definition"""

    def __init__(self, a, b, c, alpha, beta, gamma):

        self.a = a
        self.b = b
        self.c = c

    @check_sizes(points=(3, -1))
    def fractional_to_cartesian(self, points: np.ndarray):
        """Convert a list of points  from the fractional (ijk) type to cartesian (xyz)"""


    def cartesian_to_fractional(self, points: np.ndarray):
        """Convert a list of points from cartesian (xyz) to  fractional (ijk)"""
