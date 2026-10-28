"""Conversions and conventions for working with spacegroup string data."""

import re

def canonise_string(name: str):
    """Make name into a form for searching.

    Parameters
    ----------
    name : str
        Spacegroup name string to canonise.

    Returns
    -------
    str
        Lowercase string with all whitespace removed.
    """
    # remove spaces
    name = re.sub(r"\s+", "", name)

    # make lowercase
    return name.lower()

