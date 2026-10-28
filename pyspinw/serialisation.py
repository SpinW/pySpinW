"""Serialisation mixin.

See notes document 004 for details.
"""
import inspect
import os
import pathlib
from fractions import Fraction
from functools import wraps

import numpy as np
import json

from pyspinw.versions import get_version

class SPWSerialisationError(Exception):
    """Exceptions thrown by the [de]serialisation process."""


def expects_keys(keystring: str, parameter_index=0):
    """Decorator for functions that take json dicts with given keys as the `parameter_index` place parameter.

    Parameters
    ----------
    keystring : str
        Comma-separated string of required dictionary keys.
    parameter_index : int, default=0
        Positional index of the dictionary argument in the decorated function.

    Returns
    -------
    callable
        Decorator function that validates the required keys.
    """
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
                    found_keys = ",".join([key for key in json])
                    raise SPWSerialisationError(f"Expected to find key {key} in first parameter, "
                                                f"only found {found_keys}")

            return fun(*args, **kwargs)

        return wrapper
    return decorator


class SPWSerialisationContextGroup:
    """A group of objects held for serialisation.

    Parameters
    ----------
    name : str
        Name given to this context group.
    """

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
        """Generate a new ID in this context.

        Returns
        -------
        int
            Newly generated integer ID.
        """
        self._counter += 1
        return self._counter - 1

    def has(self, key):
        """Check whether this serialisation context has the object specified by `key`.

        Parameters
        ----------
        key : hashable
            Object key to look up.

        Returns
        -------
        bool
            ``True`` if `key` is registered in this context group.
        """
        return key in self._ids

    def put(self, key, serialisation_data):
        """Put a new serialised object into this context.

        Parameters
        ----------
        key : hashable
            Object key used to identify the object.
        serialisation_data : dict
            Serialised dictionary representation of the object.
        """
        if key not in self._ids:
            id = self._next_id()
            self._ids[key] = id
            self._data[id] = serialisation_data

    def reference(self, key):
        """Return a reference to an object in this group.

        Parameters
        ----------
        key : hashable
            Object key to reference.

        Returns
        -------
        dict
            Reference dictionary containing ``"reference"`` and ``"id"``.
        """
        return {"reference": self._name,
                "id": str(self._ids[key])}

    def serialise(self):
        """Get the serialisation data for this group.

        Returns
        -------
        dict
            Mapping from stringified object IDs to serialised data.
        """
        return {str(key): self._data[key] for key in self._data}



class SPWSerialisationContext:
    """Context object (see notes document 004) for serialisation.

    Attributes
    ----------
    sites : SPWSerialisationContextGroup
        Context group holding serialised lattice sites.
    """

    def __init__(self):
        self.sites = SPWSerialisationContextGroup("sites")

    def serialise(self):
        """Serialise the data in this context to JSON.

        Returns
        -------
        dict
            Dictionary of serialised context groups.
        """
        return {
            "sites": self.sites.serialise()
        }

class SPWDeserialisationRequestResponse:
    """Return type for deserialisation requests from deserialisation context.

    Contains a value which can be a deserialised object or json. Has a flag for
    whether or not it is the deserialised object.

    Parameters
    ----------
    value : object or dict
        Deserialised object instance or raw JSON dictionary.
    id : int or str
        Serialisation identifier of the requested object.
    deserialised : bool
        ``True`` if `value` has already been deserialised into an object.

    Attributes
    ----------
    value : object or dict
        Deserialised object instance or raw JSON dictionary.
    id : int or str
        Serialisation identifier of the requested object.
    deserialised : bool
        ``True`` if `value` has already been deserialised into an object.
    """

    def __init__(self, value, id: int, deserialised: bool):
        self.value = value
        self.deserialised = deserialised
        self.id = id


class SPWDeserialisationContexGroup:
    """Group of objects for deserialisation.

    Parameters
    ----------
    name : str
        Name of this context group.
    context_data : dict
        Mapping from object IDs to raw JSON representations.

    Attributes
    ----------
    objects : dict
        Collection of resolved object instances indexed by their ID.
    context_data : dict
        Mapping from object IDs to raw JSON representations.
    """

    def __init__(self, name: str, context_data: dict):
        self._name = name

        # Collection of objects that have been created, indexed by their request
        self.objects = {}
        self.context_data = context_data

    def request_by_id(self, id) -> SPWDeserialisationRequestResponse:
        """Ask for an object from this context-group using its serialisation id.

        It will return an object that contains either the deserialised object, or json that needs to
        be deserialised (along with appropriate metadata).

        Parameters
        ----------
        id : int or str
            Serialisation ID of the requested object.

        Returns
        -------
        SPWDeserialisationRequestResponse
            Response containing either the resolved instance or raw JSON data.
        """
        if id in self.objects:
            return SPWDeserialisationRequestResponse(self.objects[id], id, True)
        else:
            return SPWDeserialisationRequestResponse(self.context_data[id], id, False)

    @expects_keys("reference, id", parameter_index=1)
    def request_by_json(self, json) -> SPWDeserialisationRequestResponse:
        """Ask for an object from this context-group using a json object.

        It will return an object that contains either the deserialised object, or json that needs to
        be deserialised (along with appropriate metadata).

        Parameters
        ----------
        json : dict
            Reference dictionary containing ``"reference"`` and ``"id"`` keys.

        Returns
        -------
        SPWDeserialisationRequestResponse
            Response containing either the resolved instance or raw JSON data.

        Raises
        ------
        SPWSerialisationError
            If the reference name in `json` does not match this context group.
        """
        reference = json["reference"]
        if reference == self._name:
            return self.request_by_id(json["id"])
        else:
            raise SPWSerialisationError(f"Attempted to fetch data for '{reference}' using context for '{self.name}'")

    def put(self, id, instance):
        """Add a resolved instance to the context.

        Parameters
        ----------
        id : int or str
            Serialisation ID of the object.
        instance : object
            Resolved object instance to store.
        """
        self.objects[id] = instance

class SPWDeserialisationContext:
    """Context object (see notes document 004) for deserialisation.

    Parameters
    ----------
    context_data : dict
        Serialised context dictionary containing a ``"sites"`` mapping.

    Attributes
    ----------
    sites : SPWDeserialisationContexGroup
        Deserialisation context group for lattice sites.
    """

    @expects_keys("sites", parameter_index=1)
    def __init__(self, context_data: dict):
        self.sites = SPWDeserialisationContexGroup("sites", context_data["sites"])


def header(description: str):
    """Create header metadata for serialised files.

    Parameters
    ----------
    description : str
        User-provided description of the serialised object.

    Returns
    -------
    dict
        Metadata dictionary containing creator name, version, and description.
    """
    return {
        "created-with": "pySpinW",
        "version": get_version(),
        "description": description
    }

def save(object: "SPWSerialisable", filename: str | pathlib.Path, description: str = "", add_extension: bool = True):
    """Save a pySpinW object.

    Parameters
    ----------
    object : SPWSerialisable
        Serialisable pySpinW object to save.
    filename : str or pathlib.Path
        Destination file path.
    description : str, default=''
        Optional description stored in the file metadata.
    add_extension : bool, default=True
        If ``True``, append ``.spw`` when `filename` has no extension.
    """
    object.save(filename, description, add_extension)

def load(filename: str | pathlib.Path):
    """Load a pySpinW object from a file.

    Parameters
    ----------
    filename : str or pathlib.Path
        Path to the ``.spw`` file (with or without the ``.spw`` extension).

    Returns
    -------
    SPWSerialisable
        Deserialised pySpinW object.

    Raises
    ------
    FileNotFoundError
        If neither `filename` nor ``filename + ".spw"`` exists.
    ValueError
        If the JSON file is missing a ``"type"`` key or specifies an unknown
        object type.
    """
    # Convert to string (refactoring option would be rewriting it all in path terms)
    if isinstance(filename, pathlib.Path):
        filename = str(filename)

    with_extension = filename + ".spw"

    # Try the provided filename, if it's not there, try with extension, otherwise error
    if os.path.exists(filename) and os.path.isfile(filename):
        target_file = filename
    elif os.path.exists(with_extension) and os.path.isfile(with_extension):
        target_file = with_extension
    else:
        raise FileNotFoundError(f"Could not find file '{filename}' or '{with_extension}'")

    with open(target_file, 'r') as file:
        json_string = file.read()

    json_data = json.loads(json_string)

    try:
        serialisation_name = json_data["type"]
    except Exception:
        raise ValueError("Expected to have a 'type' key")

    from pyspinw.deserialisation import serialisation_class_lookup

    try:
        cls = serialisation_class_lookup[serialisation_name]
    except Exception:
        raise ValueError(f"Unknown object type '{serialisation_name}'")

    return cls._deserialise_from_json(json_data)


class SPWSerialisable:
    """Classes that are serialisable to SPW files should implement this interface."""

    serialisation_name = "<not-implemented>"

    def serialise(self, description: str = "") -> str:
        """Serialise an object of this type to a JSON string.

        Parameters
        ----------
        description : str, default=''
            Optional description stored in the metadata header.

        Returns
        -------
        str
            Formatted JSON string representing the object.
        """
        context = SPWSerialisationContext()
        data = {
            "meta": header(description),
            "type": self.serialisation_name,
            "object": self._serialise(context),
            "context": context.serialise()
        }
        return json.dumps(data, indent=4, sort_keys=True)

    def save(self, filename: str | pathlib.Path, description: str = "", add_extension: bool=True):
        """Save this object to a file.

        Parameters
        ----------
        filename : str or pathlib.Path
            Destination file path.
        description : str, default=''
            Optional description stored in the metadata header.
        add_extension : bool, default=True
            If ``True``, append ``.spw`` when `filename` has no extension.
        """
        filename = pathlib.Path(filename)

        if add_extension and "." not in filename.name:
            filename = filename.parent / (filename.name + ".spw")

        with open(filename, 'w') as file:
            file.write(self.serialise(description))


    @classmethod
    def deserialise(cls, json_string: str):
        """Deserialise an object of this type from a JSON string.

        Parameters
        ----------
        json_string : str
            JSON string representing the serialised object.

        Returns
        -------
        SPWSerialisable
            Deserialised instance of `cls`.
        """
        json_data = json.loads(json_string)
        return cls._deserialise_from_json(json_data)

    @classmethod
    @expects_keys("meta,type,object,context", parameter_index=1)
    def _deserialise_from_json(cls, json_data):
        """Deserialise an object of this type from a json object.

        Parameters
        ----------
        json_data : dict
            Top-level parsed JSON dictionary containing ``"meta"``, ``"type"``,
            ``"object"``, and ``"context"``.

        Returns
        -------
        SPWSerialisable
            Deserialised instance of `cls`.

        Raises
        ------
        SPWSerialisationError
            If ``json_data["type"]`` does not match ``cls.serialisation_name``.
        """
        got_type = json_data["type"]
        if cls.serialisation_name != got_type:
            raise SPWSerialisationError(f"Tried to deserialise object of kind '{got_type}' but "
                                        f"expected {cls.serialisation_name}")

        context = SPWDeserialisationContext(json_data["context"])

        return cls._deserialise(json_data["object"], context)

    def _serialise(self, context: SPWSerialisationContext) -> dict:
        """Serialise this instance into a dictionary.

        Parameters
        ----------
        context : SPWSerialisationContext
            Serialisation context for shared objects.

        Returns
        -------
        dict
            Serialised dictionary representation of this instance.

        Raises
        ------
        NotImplementedError
            If the subclass does not implement serialisation.
        """
        raise NotImplementedError("Serialisation not implemented")

    @staticmethod
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        """Deserialise an instance from a dictionary and context.

        Parameters
        ----------
        json : dict
            Serialised object dictionary.
        context : SPWDeserialisationContext
            Deserialisation context for shared objects.

        Returns
        -------
        SPWSerialisable
            Deserialised object instance.

        Raises
        ------
        NotImplementedError
            If the subclass does not implement deserialisation.
        """
        raise NotImplementedError("Deserialisation not implemented")

def rgb_serialise(rgb: tuple[float, float, float]):
    """Serialise an RGB color.

    Parameters
    ----------
    rgb : tuple of (float, float, float)
        RGB colour tuple ``(r, g, b)``.

    Returns
    -------
    dict
        Dictionary with keys ``"r"``, ``"g"``, and ``"b"``.
    """
    return {
        "r": rgb[0],
        "g": rgb[1],
        "b": rgb[2]
    }

@expects_keys("r,g,b")
def rgb_deserialise(json):
    """Deserialise an RGB color.

    Parameters
    ----------
    json : dict
        Dictionary with keys ``"r"``, ``"g"``, and ``"b"``.

    Returns
    -------
    tuple of (float, float, float)
        RGB colour tuple ``(r, g, b)``.
    """
    return (float(json["r"]),
            float(json["g"]),
            float(json["b"]))

def vec3_serialise(x,y,z):
    """Serialise a 3D vector.

    Parameters
    ----------
    x : float
        First component.
    y : float
        Second component.
    z : float
        Third component.

    Returns
    -------
    dict
        Dictionary with keys ``"x"``, ``"y"``, and ``"z"``.
    """
    return {
        "x": x,
        "y": y,
        "z": z
    }

@expects_keys("x,y,z")
def vec3_deserialise(json):
    """Deserialise a 3D vector.

    Parameters
    ----------
    json : dict
        Dictionary with keys ``"x"``, ``"y"``, and ``"z"``.

    Returns
    -------
    tuple
        3-element tuple ``(x, y, z)``.
    """
    return (json["x"], json["y"], json["z"])

def numpy_serialise(data: np.ndarray) -> dict:
    """Serialise numpy array data into a JSON-compatible dictionary.

    Parameters
    ----------
    data : numpy.ndarray
        NumPy array to serialise.

    Returns
    -------
    dict
        Dictionary containing ``"shape"``, ``"data"``, ``"imag_data"``, and
        ``"dtype"``.
    """
    shape = data.shape

    return {
        "shape": list(shape),
        "data": np.real(data).reshape(-1).tolist(),
        "imag_data": np.imag(data).reshape(-1).tolist() if np.iscomplexobj(data) else None,
        "dtype": data.dtype.str
    }


def numpy_deserialise(json: dict) -> np.ndarray:
    """Deserialise numpy array data from an SPW JSON dictionary.

    Parameters
    ----------
    json : dict
        Dictionary containing ``"shape"``, ``"data"``, ``"imag_data"``, and
        ``"dtype"``.

    Returns
    -------
    numpy.ndarray
        Reconstructed NumPy array.

    Raises
    ------
    SPWSerialisationError
        If required keys are missing, data types are invalid, or array reshaping
        fails.
    """
    try:
        shape = tuple(json["shape"])
        dtype = np.dtype(json["dtype"])

        # Don't coerce floats to ints, instead give an error
        if np.issubdtype(dtype, np.integer):
            if not np.all([isinstance(x, int) for x in json["data"]]):
                raise SPWSerialisationError("Tried to make integer numpy array from non-integer data")

        if json["imag_data"] is None:
            data = np.array(json["data"], dtype=dtype)
        else:
            data = np.array(json["data"], dtype=dtype) + 1j*np.array(json["imag_data"])

    except KeyError as ke:
        raise SPWSerialisationError(f"Failed to deserialise numpy object, bad json keys") from ke

    except ValueError as ve:
        raise SPWSerialisationError(f"Failed to deserialise values to {dtype}") from ve

    # Put in correct shape
    try:
        data = data.reshape(shape)

    except ValueError as ve:
        required_length = shape
        actual_length = len(json["data"])
        raise SPWSerialisationError(f"Failed to set shape of numpy object. "
                                    f"{shape} requires length {required_length} "
                                    f"but found {actual_length} values") from ve

    return data

def serialise_fraction(fraction: Fraction) -> dict:
    """Serialise a Fraction.

    Parameters
    ----------
    fraction : fractions.Fraction
        Fraction instance to serialise.

    Returns
    -------
    dict
        Dictionary with ``"numerator"`` and ``"denominator"`` keys.
    """
    return {
        "numerator": fraction.numerator,
        "denominator": fraction.denominator
    }

@expects_keys("numerator, denominator")
def deserialise_fraction(json: dict) -> Fraction:
    """Deserialise a Fraction.

    Parameters
    ----------
    json : dict
        Dictionary with ``"numerator"`` and ``"denominator"`` keys.

    Returns
    -------
    fractions.Fraction
        Reconstructed :class:`~fractions.Fraction`.
    """
    return Fraction(json["numerator"], json["denominator"])

def serialise_fraction_or_builtin(value: int | float | Fraction):
    """Serialise a union type of fraction or builtin number to pre-json.

    Parameters
    ----------
    value : int or float or fractions.Fraction
        Numeric value to serialise.

    Returns
    -------
    dict
        Dictionary with ``"type"`` (``"fraction"`` or ``"builtin"``) and
        ``"data"`` keys.
    """
    if isinstance(value, Fraction):
        return {
            "type": "fraction",
            "data": serialise_fraction(value)
        }
    else:
        return {
            "type": "builtin",
            "data": value
        }

@expects_keys("type, data")
def deserialise_fraction_or_builtin(json: dict) -> int | float | Fraction:
    """Deserialise the pre-json representation of a fraction or builtin.

    Parameters
    ----------
    json : dict
        Dictionary with ``"type"`` and ``"data"`` keys.

    Returns
    -------
    int or float or fractions.Fraction
        Reconstructed numeric value.

    Raises
    ------
    SPWSerialisationError
        If ``json["type"]`` is neither ``"fraction"`` nor ``"builtin"``.
    """
    match json["type"]:
        case "fraction":
            return deserialise_fraction(json["data"])
        case "builtin":
            return json["data"]
        case _:
            raise SPWSerialisationError("Expected number type to be 'fraction' or 'builtin")
