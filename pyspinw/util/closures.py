import numpy as np
from pyspinw.util.group_generators import Generator

def closure(generators: list[Generator]) -> list[Generator]:
    """ Closure of magnetic space groups"""

    output_generators = [Generator(np.eye(3), np.zeros(3), 1, name="e")]



    return output_generators