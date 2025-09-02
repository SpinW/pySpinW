""" Serialisation mixin
"""

import numpy as np

class SPWSerialisationError(Exception):
    pass

class SPWSerialisationContext:
    pass

class SPWSerialisable:

    def serialise(self, context: SPWSerialisationContext) -> dict:
        raise NotImplementedError("Serialisation not implemented")

    @staticmethod
    def deserialise(data: dict, context: SPWSerialisationContext):
        raise NotImplementedError("Deserialisation not implemented")

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