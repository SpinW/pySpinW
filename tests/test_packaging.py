""" Tests for things related to packaging """

import types

def test__all__is_present():
    """ Check that the names in __all__ in the __init__ are available the imports """
    import pyspinw

    exported = set(pyspinw.__all__)

    actual = {
        name for name in dir(pyspinw)
        if not name.startswith("_") # Ignore private
           and not isinstance(getattr(pyspinw, name), types.ModuleType) # Ignore modules
    }

    assert exported == actual, "__all__ should match imports in __init__"