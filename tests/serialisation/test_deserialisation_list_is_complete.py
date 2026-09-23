import importlib
import inspect
import pkgutil


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

def test_all_SPWSerialisable_is_represented():
    """ Gets a list of all SPWSerialisable classes, and checks there is an entry in the deserialisation list

    We need to make sure that there is an entrypoint for every `serialisation_name`
    """

    classes = get_classes("pyspinw")

    print(classes)