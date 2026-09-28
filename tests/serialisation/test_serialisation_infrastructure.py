import importlib
import inspect
import pkgutil

from pyspinw.serialisation import SPWSerialisable
from pyspinw.deserialisation import serialisation_class_lookup


def get_classes(package_name):
    package = importlib.import_module(package_name)

    classes = []

    # Include the package itself (__init__.py)
    modules = [package]

    # Include all submodules/subpackages
    if hasattr(package, "__path__"):
        for module_info in pkgutil.walk_packages(
            package.__path__,
            package.__name__ + "."
        ):
            modules.append(importlib.import_module(module_info.name))

    for module in modules:
        for name, obj in inspect.getmembers(module, inspect.isclass):
            classes.append((module.__name__, name, obj))

    return classes

def test_all_SPWSerialisable_is_implemented():
    """ Check that no class other than the SPWSerialisable mixin have a "<not-implemented>" serialisation name """

    # First we get all the classes
    classes = get_classes("pyspinw")

    # Then filter for the serialisable classes
    serialisable_classes = [cls for mod, name, cls in classes
                            if issubclass(cls, SPWSerialisable)]

    for cls in serialisable_classes:
        # Skip the base class
        if "SPWSerialisable" in cls.__name__:
            continue

        assert cls.serialisation_name != "<not-implemented>", f"{cls.__name__} should have a serialisation name"

def test_all_serialisation_names_follow_conventions():
    """ Check that serialisation names follow our conventions """

    # First we get all the classes
    classes = get_classes("pyspinw")

    # Then filter for the serialisable classes
    serialisable_classes = [cls for mod, name, cls in classes
                            if issubclass(cls, SPWSerialisable)]

    for cls in serialisable_classes:
        # Skip the base class
        if "SPWSerialisable" in cls.__name__:
            continue

        assert cls.serialisation_name == cls.serialisation_name.lower(), \
            f"`serialisation_name` should be lowercase only, got '{cls.serialisation_name}'"


        assert " " not in cls.serialisation_name, \
            f"`serialisation_name` should not have spaces, got '{cls.serialisation_name}'"

        for forbidden in "<>-:/\[]":
            assert forbidden not in cls.serialisation_name, \
                f"`serialisation_name` should not have '{forbidden}', got '{cls.serialisation_name}'"


def test_all_SPWSerialisable_is_represented():
    """ Gets a list of all SPWSerialisable classes, and checks there is an entry in the deserialisation list

    We need to make sure that there is an entrypoint for every `serialisation_name`
    """


    # First we get all the classes
    classes = get_classes("pyspinw")

    # Then filter for the serialisable classes
    serialisable_classes = [cls
                            for mod, name, cls in classes
                            if issubclass(cls, SPWSerialisable)]

    # Get a set of the serialisation_names
    names = set(cls.serialisation_name for cls in serialisable_classes if "SPWSerialisable" not in cls.__name__)

    # Check the deserialisation list
    for name in names:
        assert name in serialisation_class_lookup, f"There should be an entry in deserialisation.py for `{name}`"