from dataclasses import dataclass


@dataclass
class Tolerances:
    SAME_SITE_ABS_TOL = 1e-10
    VECTOR_TOL = 1e-10
    COUPLING_ORDER_THRESHOLD = 1e-7
    IS_INTEGER_TOL = 1e-6
    IS_ZERO_TOL = 1e-10

tolerances = Tolerances()