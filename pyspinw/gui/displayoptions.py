"""Dataclass for display/rendering options."""

import json
from dataclasses import dataclass, asdict


@dataclass
class DisplayOptions:
    """Options for how things should display in the viewer.

    Parameters
    ----------
    show_sites : bool, optional
        Whether to render atomic/magnetic sites. Default is ``True``.
    show_exchanges : bool, optional
        Whether to render exchange interactions. Default is ``True``.
    show_anisotropies : bool, optional
        Whether to render single-ion anisotropies. Default is ``True``.
    show_unit_cell : bool, optional
        Whether to render unit cell wireframes. Default is ``True``.
    show_supercell : bool, optional
        Whether to render the supercell wireframe. Default is ``True``.
    prettify : bool, optional
        Whether to render periodic boundary images on cell faces, edges, and
        corners. Default is ``True``.
    show_nonmagnetic_atoms : bool, optional
        Whether to render non-magnetic sites. Default is ``True``.
    use_atomic_radii : bool, optional
        Whether to scale site spheres by atomic radius. Default is ``True``.
    show_spins_not_atoms : bool, optional
        Whether to render magnetic sites as spin arrows rather than spheres.
        Default is ``True``.
    atom_spin_scaling : float, optional
        Scaling factor for site spheres and spin arrows. Default is ``0.35``.
    exchange_scaling : float, optional
        Scaling factor for exchange bond tubes. Default is ``0.20``.
    anisotropy_scaling : float, optional
        Scaling factor for anisotropy ellipsoids. Default is ``0.35``.
    show_cartesian_axes : bool, optional
        Whether to display the Cartesian orientation axes. Default is ``True``.
    show_lattice_axes : bool, optional
        Whether to display the crystal lattice orientation axes. Default is
        ``False``.
    orthogonal_lattice_axes : bool, optional
        Whether to orthogonalise the displayed lattice axes. Default is
        ``False``.
    perspective : bool, optional
        Whether to use perspective projection instead of orthographic. Default
        is ``True``.
    background_color : tuple of (float, float, float), optional
        RGB background color. Default is ``(0.05, 0.05, 0.08)``.
    default_exchange_color : tuple of (float, float, float), optional
        Default RGB color for exchange bonds. Default is ``(0.2, 0.4, 0.8)``.
    default_site_color : tuple of (float, float, float), optional
        Default RGB color for sites. Default is ``(0.7, 0.8, 0.6)``.
    selected_color : tuple of (float, float, float), optional
        RGB highlight color for selected elements. Default is
        ``(1.0, 0.6, 0.1)``.
    hover_color : tuple of (float, float, float), optional
        RGB highlight color for hovered elements. Default is
        ``(1.0, 1.0, 1.0)``.
    selected_hover_color : tuple of (float, float, float), optional
        RGB highlight color for elements that are both selected and hovered.
        Default is ``(1.0, 0.8, 0.1)``.
    anisotropy_color : tuple of (float, float, float), optional
        RGB color for anisotropy ellipsoids. Default is ``(0.2, 0.5, 1.0)``.
    """

    show_sites: bool = True
    show_exchanges: bool = True
    show_anisotropies: bool = True
    show_unit_cell: bool = True
    show_supercell: bool = True

    prettify: bool = True

    show_nonmagnetic_atoms: bool = True
    use_atomic_radii: bool = True
    show_spins_not_atoms: bool = True

    atom_spin_scaling: float = 0.35
    exchange_scaling: float = 0.20
    anisotropy_scaling: float = 0.35

    show_cartesian_axes: bool = True
    show_lattice_axes: bool = False
    orthogonal_lattice_axes: bool = False

    perspective: bool = True

    background_color: tuple[float, float, float] = 0.05, 0.05, 0.08
    default_exchange_color: tuple[float, float, float] = 0.2, 0.4, 0.8
    default_site_color: tuple[float, float, float] = 0.7, 0.8, 0.6

    selected_color: tuple[float, float, float] = 1.0, 0.6, 0.1
    hover_color: tuple[float, float, float] = 1.0, 1.0, 1.0
    selected_hover_color: tuple[float, float, float] = 1.0, 0.8, 0.1

    anisotropy_color: tuple[float, float, float] = 0.2, 0.5, 1.0

    def serialise(self) -> str:
        """Serialise settings object to a string.

        Returns
        -------
        str
            JSON-serialised string representation of the display options.
        """
        out = asdict(self)

        return json.dumps(out)

    @staticmethod
    def deserialise(serialised: str):
        """Deserialise a setting object from a string.

        Parameters
        ----------
        serialised : str
            JSON-serialised string representation of the display options.

        Returns
        -------
        DisplayOptions
            Deserialised ``DisplayOptions`` instance, or default options if
            deserialisation fails.
        """
        try:
            data = json.loads(serialised)

            return DisplayOptions(**data)

        except Exception: # Anything at all goes wrong, just return default

            return DisplayOptions()

    @staticmethod
    def LightMode() -> 'DisplayOptions':
        """Returns a set of viewer ``DisplayOptions`` with a light background"""
        display_options = DisplayOptions()
        display_options.background_color = 0.95, 0.95, 0.92
        display_options.default_site_color = 0.5, 0.6, 0.4
        display_options.hover_color = 0.0, 0.0, 0.0
        display_options.selected_color = 0.8, 0.5, 0.1
        return display_options

    @staticmethod
    def DarkMode() -> 'DisplayOptions':
        """Returns a set of viwer ``DisplayOptions`` with a dark background"""
        return DisplayOptions()
