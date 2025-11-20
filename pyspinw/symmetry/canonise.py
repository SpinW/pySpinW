""" Conversions and conventions for working with spacegroup string data"""

import re

def canonise_string(name: str):
    """ Make name into a form for searching """
    # remove spaces
    name = re.sub(r"\s+", "", name)

    # make lowercase
    return name.lower()

