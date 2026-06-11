""" Dataclass for display/rendering options"""

import json
from dataclasses import dataclass, asdict


@dataclass
class DisplayOptions:
    """ Options for how things should display """

    show_sites: bool = True
    show_exchanges: bool = True
    show_anisotropies: bool = True
    show_unit_cell: bool = False
    show_supercell: bool = False

    prettify: bool = True

    show_nonmagnetic_atoms: bool = True
    use_atomic_radii: bool = True
    show_atoms_not_spins: bool = False

    atom_spin_scaling: float = 0.35
    exchange_scaling: float = 0.20

    show_cartesian_axes: bool = True
    show_lattice_axes: bool = False
    orthogonal_lattice_axes: bool = False

    background_color: tuple[float, float, float] = 0.05, 0.05, 0.08
    default_exchange_color: tuple[float, float, float] = 0.2, 0.4, 0.8
    default_site_color: tuple[float, float, float] = 0.7, 0.8, 0.6

    selected_color: tuple[float, float, float] = 1.0, 0.6, 0.1
    hover_color: tuple[float, float, float] = 1.0, 1.0, 1.0
    selected_hover_color: tuple[float, float, float] = 1.0, 0.8, 0.1

    def serialise(self) -> str:
        """ Serialise settings object to a string"""
        out = asdict(self)

        return json.dumps(out)

    @staticmethod
    def deserialise(serialised: str):
        """ Deserialise a setting object from a string"""
        try:
            data = json.loads(serialised)

            return DisplayOptions(**data)

        except Exception: # Anything at all goes wrong, just return default

            return DisplayOptions()
