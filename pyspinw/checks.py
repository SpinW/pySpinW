"""Tools for checking input shapes"""

import functools
import inspect

from typing import Callable
from collections import defaultdict

import numpy as np


class DimensionalityError(ValueError):
    """The dimensions of a numpy array don't match the specification"""


def check_sizes(force_numpy: bool = False, allow_nones: bool = False, **kwargs):
    """Decorator to check the dimensionality of a given vector

    Example usage:
        check that a function has an n-by-3 vector input called `a`, and
        a 2-by-n-by-whatever input called 'b'

        @dimensionality_check(a=('n',3), b=(2,'n',-1))
        def my_function(a, b):
            ...

        Limitations: can only check dimensions are the same, can't check things like whether sizes are related by
                     some arbitrary formula.

    :param force_numpy: default=False, Convert the named arrays into numpy form if they are not already
    :param allow_nones: default=False, Allow None to be given for values

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

    for keyword, constraint in kwargs.items():

        if not isinstance(constraint, tuple):
            raise TypeError("Size constrains should be specified by tuple of int/str. "
                            f"Got type {type(constraint)} for '{keyword}'")

        sizes.append((keyword, len(constraint)))

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
        variable_names = fun.__code__.co_varnames

        sig = inspect.signature(fun)
        defaults = {
            name: param.default
            for name, param in sig.parameters.items()
            if param.default is not inspect._empty
        }

        @functools.wraps(fun)
        def wrapper(*args, **kwargs):

            # shove all the arguments in a dictionary keyed by the variable names
            all_args = defaults.copy()
            for name, arg in zip(variable_names, args):
                all_args[name] = arg
            all_args.update(**kwargs)

            # Allow skipping of size checks
            if "skip_size_check" in kwargs and all_args["skip_size_check"]:
                del all_args["skip_size_check"]
                return fun(**all_args)

            # Check the sizes, and the type while were at it
            #
            #  Note: THIS POTENTIALLY UPDATES all_args
            #
            for name, size in sizes:
                if allow_nones and all_args[name] is None:
                    continue

                if name not in all_args:
                    raise ValueError(f"The required numpy array '{name}' was not given")

                if not isinstance(all_args[name], np.ndarray):
                    if force_numpy:
                        all_args[name] = np.array(all_args[name])
                    else:
                        raise TypeError(f"Argument '{name}' is not a numpy array, but is {type(all_args[name])}")

                if not len(all_args[name].shape) == size:
                    raise DimensionalityError(f"Expected '{name}' to be a {size}D tensor, "
                                              f"but it is {len(all_args[name].shape)}D")

            # Check all the constant values
            for name, dimension, size in constants:
                data = all_args[name]

                if allow_nones and data is None:
                    continue

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
                    if allow_nones and data is None:
                        continue

                    # Again, the validity of this should be checked in len(shape) loop
                    if size is None:
                        size = data.shape[dimension]
                    else:
                        if not size == data.shape[dimension]:
                            message_details = ", ".join([f"{name}:{index}" for name, index in equalities[symbol]])
                            raise DimensionalityError(f"Tensor dimensions specified by '{symbol}' do not all match "
                                                      f"({message_details})")

            # if all these have passed then we're peachy

            # Call function with potentially updated arguments, by keyword, even though we used args before
            return fun(**all_args)

        return wrapper

    return decorator
