""" Serialisation mixin
"""
import inspect
from functools import wraps

import numpy as np

class SPWSerialisationError(Exception):
    pass


class SPWSerialisationContextGroup:
    def __init__(self, name: str):
        # Name given to this context
        self._name = name

        # For generating keys
        self._counter = 0

    def _next_id(self):
        """ Generate a new ID in this context """
        self._counter += 1
        return self._counter - 1


class SPWSerialisationContext:
    def __init__(self):
        self.sites = SPWSerialisationContextGroup("site")

class SPWDeserialisationContexGroup:
    def __init__(self, name: str):
        self._name = name

class SPWDeserialisationContext:
    def __init__(self):
        self.sites = SPWDeserialisationContexGroup("site")

class SPWSerialisable:
    """ Classes that are serialisable to SPW files should use implement this interface """

    def serialise(self):
        pass

    @staticmethod
    def deserialise(json):
        pass

    def _serialise(self, context: SPWSerialisationContext) -> dict:
        raise NotImplementedError("Serialisation not implemented")

    @staticmethod
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        raise NotImplementedError("Deserialisation not implemented")


def expects_keys(keystring: str):
    """ Decorator for functions that take json dicts as first argument with given keys"""
    keys = [x.strip() for x in keystring.split(",")]

    def decorator(fun):

        # We want to get the first parameter, regardless of how it is called, need to use inspection
        sig = inspect.signature(fun)
        first_param = next(iter(sig.parameters))

        @wraps(fun)
        def wrapper(*args, **kwargs):

            # Get our arguments by name, even if they're positional
            bound = sig.bind(*args, **kwargs)
            bound.apply_defaults()
            json = bound.arguments[first_param]

            if not isinstance(json, dict):
                raise SPWSerialisationError(f"Expected a dictionary as first argument, got {json}")

            for key in keys:
                if key not in json:
                    raise SPWSerialisationError(f"Expected to find key {key} in first parameter")

            return fun(*args, **kwargs)

        return wrapper
    return decorator

def numpy_serialise(data: np.ndarray) -> dict:
    shape = data.shape

    return {
        "shape": list(shape),
        "data": data.reshape(-1).tolist(),
        "dtype": data.dtype.str
    }


def numpy_deserialise(json: dict) -> np.ndarray:
    try:
        shape = tuple(json["shape"])
        dtype = np.dtype(json["dtype"])

        # Don't coerce floats to ints, instead give an error
        if np.issubdtype(dtype, np.integer):
            if not np.all([isinstance(x, int) for x in json["data"]]):
                raise SPWSerialisationError("Tried to make integer numpy array from non-integer data")

        data = np.array(json["data"], dtype=dtype)

    except KeyError as ke:
        raise SPWSerialisationError(f"Failed to deserialise numpy object, bad json keys") from ke

    except ValueError as ve:
        raise SPWSerialisationError(f"Failed to deserialise values to {dtype}") from ve

    # Put in correct shape
    try:
        data = data.reshape(shape)

    except ValueError as ve:
        required_length = shape
        raise SPWSerialisationError(f"Failed to set shape of numpy object. "
                                    f"{shape} requires length {required_length} "
                                    f"but found {len(json["data"])} values") from ve

    return data
