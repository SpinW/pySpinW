"""Data for working with tolerances"""

from dataclasses import dataclass


@dataclass
class Tolerances:
    """Holds different tolerance values for numerical methods """

    SAME_SITE_ABS_TOL = 1e-10
    VECTOR_TOL = 1e-10
    EXCHANGE_ORDER_THRESHOLD = 1e-7
    IS_INTEGER_TOL = 1e-6
    IS_ZERO_TOL = 1e-10
    BOND_TOL = 1e-6        # Distance tolerance to consider same bond
    SPECIALISE_ROUNDING_EXPONENT = 14
    IMAG_MODE_TOL = 1e-2   # Value of imaginary energy (meV) below which to ignore

tolerances = Tolerances()
