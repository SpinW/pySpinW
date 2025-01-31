import numpy as np

Generator = tuple[np.ndarray, np.ndarray, float]

def compose_generators(a: Generator, b: Generator) -> Generator:
    """ Calculate the composition x -> a(b(x) """

    #
    # for the vector part
    # a(x) = Ma x + Va
    # b(x) = Mb x + Vb
    #
    # a(b(x)) = Ma(Mb + Vb) + Va
    #         = (Ma Mb) + (Ma Vb + Va)
    #

    return a[0] @ b[0], a[0] @ b[1] + a[1], a[2] * b[2]


def order_generators(list[Generator]) -> list[Generator]:
    """ Put generators in order """

    # Convert to big list, 9 + 3 + 1 = 13 across
    

def closure(generators: list[Generator]) -> list[Generator]:
    """ Closure of magnetic space groups"""


    pass