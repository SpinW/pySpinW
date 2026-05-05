""" Metadata for sites """

import re

from ase.data.colors import jmol_colors
from ase.data import atomic_numbers, chemical_symbols, covalent_radii, missing

from pyspinw.serialisation import SPWSerialisable, SPWSerialisationContext, SPWDeserialisationContext, rgb_serialise, \
    expects_keys, rgb_deserialise


def _build_regex():
    """ Makes the regex for matching elements"""
    # Sort so 2-letter symbols come before 1-letter ones
    symbols_sorted = sorted(chemical_symbols[1:], key=lambda x: (-len(x), x))

    # Build regex
    pattern = r"^(?:{})".format("|".join(symbols_sorted))
    return re.compile(pattern)


_chemical_pattern = _build_regex()

class SiteMetadata(SPWSerialisable):
    """ Metadata about the element in the site """

    def __init__(self,
                 element: str | None = None,
                 radius: float | None = None,
                 color: tuple[float, float, float] | None=None):

        # Need to ignore entry 0 as this is "X"
        if element is not None and element not in chemical_symbols[1:]:
            raise ValueError(f"'{element}' is not a valid element")

        self.element = element
        self.radius = radius
        self.color = color

    @staticmethod
    def metadata_from_name(name: str | None):
        """ Extract metadata based on name"""
        if name is None:
            return SiteMetadata()

        # See if the name starts with a known element, the first in the ASE list is X, ignore that
        match_result = _chemical_pattern.match(name)

        if match_result:
            element = match_result.group()
            z = atomic_numbers[element]

            color = jmol_colors[z]
            radius = covalent_radii[z]
            radius = 1.0 if radius == missing else radius

            return SiteMetadata(element, radius, color)


        else:
            return SiteMetadata(None, None, None)


    def _serialise(self, context: SPWSerialisationContext) -> dict:
        return {
            "element": self.element,
            "radius": self.radius,
            "color": None if self.color is None else rgb_serialise(self.color)
        }

    @staticmethod
    @expects_keys("element, radius, color")
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        element = json["element"]
        radius = json["radius"]
        color = json["color"]
        color = color if color is None else rgb_deserialise(color)

        return SiteMetadata(element, radius, color)
