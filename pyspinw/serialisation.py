""" Serialisation mixin
"""
import inspect
from functools import wraps

import numpy as np

class SPWSerialisationError(Exception):
    pass


def expects_keys(keystring: str, parameter_index=0):
    """ Decorator for functions that take json dicts as first (or other if parameter_index is specified)
    argument with given keys"""
    keys = [x.strip() for x in keystring.split(",")]

    def decorator(fun):

        # We want to get the first parameter, regardless of how it is called, need to use inspection
        sig = inspect.signature(fun)
        params = list(sig.parameters)

        if len(params) > parameter_index:
            parameter = params[parameter_index]
        else:
            raise RuntimeError("expects_keys decorator cannot be applied with `parameter_index` out of range")

        @wraps(fun)
        def wrapper(*args, **kwargs):

            # Get our arguments by name, even if they're positional
            bound = sig.bind(*args, **kwargs)
            bound.apply_defaults()
            json = bound.arguments[parameter]

            if not isinstance(json, dict):
                raise SPWSerialisationError(f"Expected a dictionary as first argument, got {json}")

            for key in keys:
                if key not in json:
                    raise SPWSerialisationError(f"Expected to find key {key} in first parameter")

            return fun(*args, **kwargs)

        return wrapper
    return decorator


class SPWSerialisationContextGroup:
    """ A group of objects held for serialisation"""

    def __init__(self, name: str):
        # Name given to this context
        self._name = name

        # For generating keys
        self._counter = 0

        # ID lookup (object key to id)
        self._ids = {}

        # Data lookup (id to json data)
        self._data = {}

    def _next_id(self):
        """ Generate a new ID in this context """
        self._counter += 1
        return self._counter - 1

    def put(self, key, serialisation_data):
        if key not in self._ids:
            id = self._next_id()
            self._ids[key] = id
            self.data[id] = serialisation_data

    def reference(self, key):
        """ Return a reference to an object in this group"""
        return {"reference": self.name,
                "id": self._ids[key]}

    def serialise(self):
        """ Get the serialisation data for this group """
        return {str(key): self._data[key] for key in self._data}



class SPWSerialisationContext:
    def __init__(self):
        self.sites = SPWSerialisationContextGroup("sites")

    def serialise(self):
        return {
            "sites": self.sites.serialise()
        }

class SPWDeserialisationRequestResponse:
    """ Return type for deserialisation requests from deserialisation context

    Contains a value which can be a deserialised object or json. Has a flag for
    is whether or not it is the deserialised object
    """
    def __init__(self, value, deserialised: bool):
        self.value = value
        self.deserialised = deserialised


class SPWDeserialisationContexGroup:
    """ Group of objects for deserialisation """
    def __init__(self, name: str, context_data: dict):
        self._name = name

        # Collection of objects that have been created, indexed by their request
        self.objects = {}
        self.context_data = context_data

    def request_by_id(self, id) -> SPWDeserialisationRequestResponse:
        if id in self.objects:
            return SPWDeserialisationRequestResponse(self.objects[id], True)
        else:
            return SPWDeserialisationRequestResponse(self.context_data[id], False)

    @expects_keys("reference, id", parameter_index=1)
    def request_by_json(self, json) -> SPWDeserialisationRequestResponse:
        reference = json["reference"]
        if reference == self._name:
            return self.request_by_id(json["id"])
        else:
            raise SPWSerialisationError(f"Attempted to fetch data for '{reference}' using context for '{self.name}'")

    def put(self, id, instance):
        """ Add a resolved instance to the context """
        self.objects[id] = instance

class SPWDeserialisationContext:
    """ Carries information needed to deserialise objects"""

    @expects_keys("sites", parameter_index=1)
    def __init__(self, context_data: dict):
        self.sites = SPWDeserialisationContexGroup("sites", context_data["sites"])


class SPWSerialisable:
    """ Classes that are serialisable to SPW files should use implement this interface """

    serialisation_name = "<not-implemented>"

    def serialise(self):
        context = SPWSerialisationContext()
        return {
            "type": self.serialisation_name,
            "object": self._serialise(context),
            "context": context.serialise()
        }

    @classmethod
    @expects_keys("type, object, context", parameter_index=1)
    def deserialise(cls, json):
        got_type = json["type"]
        if cls.serialisation_name != got_type:
            raise SPWSerialisationError(f"Tried to deserialise object of kind '{got_type}' but "
                                        f"expected {cls.serialisation_name}")

        context = SPWDeserialisationContext(json["context"])

        cls._deserialise(json["object"], context)

    def _serialise(self, context: SPWSerialisationContext) -> dict:
        raise NotImplementedError("Serialisation not implemented")

    @staticmethod
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        raise NotImplementedError("Deserialisation not implemented")


def numpy_serialise(data: np.ndarray) -> dict:
    """ Serialisation for numpy data in the SPW json"""
    shape = data.shape

    return {
        "shape": list(shape),
        "data": data.reshape(-1).tolist(),
        "dtype": data.dtype.str
    }


def numpy_deserialise(json: dict) -> np.ndarray:
    """ Deserialisation method numpy data in SPW json"""
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
