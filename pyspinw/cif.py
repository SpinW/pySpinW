""" CIF File Loading """
import numpy as np
from CifFile import ReadCif

from pyspinw import LatticeSite, TiledSupercell, Structure
from pyspinw.sitemeta import SiteMetadata
from pyspinw.symmetry.group import NoSuchGroup
from pyspinw.symmetry.supercell import Supercell
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.interface import spacegroup


def parse_float_entry(s: str):
    """ Parse a float entry, they can have brackets for values that go beyond the specified precision"""

    s = s.replace("(", "")
    s = s.replace(")", "")

    return float(s)

def load_cif(filename: str, supercell: Supercell = TiledSupercell(), entry_index: int=0):
    """ Load a CIF file in as a structure """

    # Get the right bit of data

    cif = ReadCif(filename)

    names = cif.keys()

    try:
        name = names[entry_index]
    except IndexError:
        raise ValueError(f"Entry index out of range (file has {len(names)} entries)")

    data = cif[name]

    # Get the spacegroup and lattice parameters

    spacegroup_name = data["_symmetry_space_group_name_H-M"]

    sg = spacegroup(spacegroup_name)


    a = parse_float_entry(data["_cell_length_a"])
    b = parse_float_entry(data["_cell_length_b"])
    c = parse_float_entry(data["_cell_length_c"])

    alpha = parse_float_entry(data["_cell_angle_alpha"])
    beta = parse_float_entry(data["_cell_angle_beta"])
    gamma = parse_float_entry(data["_cell_angle_gamma"])

    cell = UnitCell(a,b,c,alpha=alpha,beta=beta,gamma=gamma)

    # Atom radii if given

    radius_lookup = None

    if "_atom_type_radius_bond" and "_atom_type_symbol" in data:
        radius_lookup = {atom: parse_float_entry(radius)
         for atom, radius in zip(data["_atom_type_symbol"], data["_atom_type_radius_bond"])}

    sites = []

    # Create the lattice sites
    for label, atom, x_string, y_string, z_string in zip(
            data["_atom_site_label"],
            data["_atom_site_type_symbol"],
            data["_atom_site_fract_x"],
            data["_atom_site_fract_y"],
            data["_atom_site_fract_z"]):

        metadata = SiteMetadata.metadata_from_name(label)
        metadata.element = atom

        if radius_lookup is not None:
            metadata.radius = radius_lookup[atom]

        x = parse_float_entry(x_string)
        y = parse_float_entry(y_string)
        z = parse_float_entry(z_string)

        # We need to set sensible defaults according to the supercell
        supercell_spins = np.zeros((supercell.n_components(), 3), dtype=float)

        site = LatticeSite(x,y,z, supercell_spins=supercell_spins, name=label, metadata=metadata)

        sites.append(site)

    return Structure(sites, unit_cell = cell, spacegroup=sg, supercell=supercell)

if __name__ == "__main__":
    structure = load_cif("../example_structures/1100231.cif")

    print(structure.site_by_name("Si1").metadata)

    from pyspinw import view
    view(structure)
