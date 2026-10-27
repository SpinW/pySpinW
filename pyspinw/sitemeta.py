"""Metadata for sites."""

import re

from ase.data.colors import jmol_colors
from ase.data import atomic_numbers, chemical_symbols, covalent_radii, missing

from pyspinw.serialisation import SPWSerialisable, SPWSerialisationContext, SPWDeserialisationContext, rgb_serialise, \
    expects_keys, rgb_deserialise


def _build_regex():
    """Make the regex for matching elements.

    Returns
    -------
    re.Pattern
        Compiled regular expression matching chemical element symbols at the
        start of a string.
    """
    # Sort so 2-letter symbols come before 1-letter ones
    symbols_sorted = sorted(chemical_symbols[1:], key=lambda x: (-len(x), x))

    # Build regex
    pattern = r"^\s*({})(?:[^a-zA-Z]+.*)?$".format("|".join(symbols_sorted))
    return re.compile(pattern)


_chemical_pattern = _build_regex()

class SiteMetadata(SPWSerialisable):
    """Metadata about the element in the site.

    Parameters
    ----------
    element : str or None, optional
        Chemical element symbol (e.g. ``"Fe"``). Default is ``None``.
    radius : float or None, optional
        Atomic radius used for rendering. Default is ``None``.
    color : tuple of (float, float, float) or None, optional
        RGB colour tuple used for rendering. Default is ``None``.

    Attributes
    ----------
    element : str or None
        Chemical element symbol.
    radius : float or None
        Atomic radius used for rendering.
    color : tuple of (float, float, float) or None
        RGB colour tuple used for rendering.
    """

    serialisation_name = "site_metadata"

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

    def __repr__(self):
        return f"SiteMetadata(element={self.element}, radius={self.radius}, color={self.color})"

    def copy(self):
        """Make a copy of this metadata.

        Returns
        -------
        SiteMetadata
            Copy of this ``SiteMetadata`` instance.
        """
        return SiteMetadata(self.element, self.radius, self.color)

    @staticmethod
    def extract_element(string: str):
        """Get a chemical element from a string.

        Parameters
        ----------
        string : str
            Input string (such as a site name) from which to extract a
            chemical symbol.

        Returns
        -------
        str or None
            Matched chemical element symbol, or ``None`` if no valid element
            prefix is found.
        """
        match_result = _chemical_pattern.match(string)

        if match_result:
            return match_result.group(1)

        else:
            return None

    @staticmethod
    def metadata_from_name(name: str | None):
        """Extract metadata based on name.

        Parameters
        ----------
        name : str or None
            Site name from which to infer element, radius, and colour.

        Returns
        -------
        SiteMetadata
            Inferred metadata for the given site name.
        """
        if name is None:
            return SiteMetadata()

        # See if the name starts with a known element, the first in the ASE list is X, ignore that
        element = SiteMetadata.extract_element(name)

        if element is not None:
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
