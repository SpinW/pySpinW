from dataclasses import dataclass


@dataclass
class Tolerances:
    SAME_SITE_ABS_TOL = 1e-10
    VECTOR_TOL = 1e-10

tolerances = Tolerances()