import numpy as np
from typing import Callable, Any
from collections import defaultdict

class DimensionalityError(ValueError):
    """ The dimensions of a numpy array don't match the specification """

def dimensionality_check(**kwargs):
    """ Decorator to check the dimensionality of a given vector

    Example usage:
        check that a function has an n-by-3 vector input called `a`, and
        a 2-by-n-by-whatever input called 'b'

        @dimensionality_check(a=('n',3), b=(2,'n',-1))
        def my_function(a, b):
            ...

        Limitations: can only check dimensions are the same, can't check things like whether sizes are related by
                     some arbitrary formula.

    """

    # If assertions are not enabled, give a decorator that just return the function as is
    if not __debug__:
        def identity(fun):
            return fun

        return identity

    #
    # We're not in debug mode, so do the work to construct the decorator that does the check
    #

    # lookup of variable name and indices that have to be the same, indexed by symbol
    equalities = defaultdict(list[tuple[str, int]])

    constants: list[tuple[str, int, int]] = [] # List of variable names and indices that are fixed
    sizes: list[tuple[str, int]] = [] # list of len(shape)

    for keyword in kwargs:
        constraint = kwargs[keyword]

        sizes.append((keyword, len(sizes)))

        for dimension, size in enumerate(constraint):

            if isinstance(size, int):
                if size >= 0:
                    constants.append((keyword, dimension, size))
                elif size < -1:
                    raise ValueError("Prescribed array size must be positive (or -1 for no constraint). "
                                     f"Got {keyword}.shape[{dimension}] = {size}")

            elif isinstance(size, str):
                equalities[size].append((keyword, dimension))

            else:
                raise TypeError("Constraints must be int or str, got "
                                f"{keyword}.shape[{dimension}] constraint of type {type(size)}")

    def decorator(fun: Callable) -> Callable:

        # grab the argument names
        variable_names = fun.__code__.co_argnames
        def wrapper(*args, **kwargs):

            # shove all the arguments in a dictionary keyed by the variable names
            all_args = {name: arg for name, arg in zip(variable_names, args)}
            all_args.update(**kwargs)

            # Check the sizes, and the type while were at it
            for name, size in sizes:

                if name not in all_args:
                    raise ValueError(f"The numpy array required ('{name}') was not given")

                data = all_args[name]

                if isinstance(data, np.ndarray):
                    if not len(data.shape) == size:
                        raise DimensionalityError(f"Expected '{name}' to be a {size}D tensor, "
                                                  f"but it is {len(data.shape)}D")
                else:
                    raise TypeError(f"Argument '{name}' is not a numpy array, but is {type(data)}")

            # Check all the constant values
            for name, dimension, size in constants:
                data = all_args[name]
                # We have already checked that it's an array and for correct len(shape), so resolving
                # this should not throw

                if not data.shape[dimension] == size:
                    raise DimensionalityError(f"Dimension {dimension} of '{name}' is not {size}. "
                                              f"Got {data.shape[dimension]}")

            # Check all the equalities
            for symbol in equalities:
                size = None
                for name, dimension in equalities[symbol]:
                    data = all_args[name]
                    # Again, the validity of this should be checked in len(shape) loop
                    if size is None:
                        size = data.shape[dimension]
                    else:
                        if not size == data.shape[dimension]:
                            message_details = ", ".join(["%s:%i" % tup for tup in equalities[symbol]])
                            raise DimensionalityError(f"Tensor dimensions specified by '{symbol}' do not all match "
                                                      f"({message_details})")

            # if all these have passed then we're peachy

            return fun(*args, **kwargs)

        return wrapper

    return decorator
