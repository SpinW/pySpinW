"""Data for working with tolerances."""

from dataclasses import dataclass


@dataclass
class Tolerances:
    """Holds different tolerance values for numerical methods.

    Attributes
    ----------
    SAME_SITE_ABS_TOL : float
        Absolute tolerance for considering two site positions identical.
    VECTOR_TOL : float
        Tolerance for vector component equality comparisons.
    EXCHANGE_ORDER_THRESHOLD : float
        Distance difference threshold for grouping exchanges into shells.
    IS_INTEGER_TOL : float
        Tolerance for checking whether a floating-point value is an integer.
    IS_ZERO_TOL : float
        Tolerance for checking whether a value or vector magnitude is zero.
    BOND_TOL : float
        Distance tolerance to consider two bonds as belonging to the same shell.
    SPECIALISE_ROUNDING_EXPONENT : int
        Decimal rounding exponent used when specialising exchanges and anisotropies.
    IMAG_MODE_TOL : float
        Value of imaginary energy (meV) below which to ignore imaginary modes.
    """

    SAME_SITE_ABS_TOL = 1e-10
    VECTOR_TOL = 1e-10
    EXCHANGE_ORDER_THRESHOLD = 1e-7
    IS_INTEGER_TOL = 1e-6
    IS_ZERO_TOL = 1e-10
    BOND_TOL = 1e-6        # Distance tolerance to consider same bond
    SPECIALISE_ROUNDING_EXPONENT = 14
    IMAG_MODE_TOL = 1e-2   # Value of imaginary energy (meV) below which to ignore

tolerances = Tolerances()
