""" CIF File Loading """
import numpy as np
from CifFile import ReadCif
import pathlib

from CifFile.CifFile_module import ReadCifWithErrors

from pyspinw import LatticeSite, TiledSupercell, Structure
from pyspinw.sitemeta import SiteMetadata
from pyspinw.symmetry.group import NoSuchGroup
from pyspinw.symmetry.magnetic_group import MagneticSpaceGroup
from pyspinw.symmetry.operations import MagneticOperation
from pyspinw.symmetry.supercell import Supercell
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.interface import spacegroup




def parse_float_entry(s: str):
    """ Parse a float entry, they can have brackets for values that go beyond the specified precision"""
    s = s.replace("(", "")
    s = s.replace(")", "")

    return float(s)

class LoadError(Exception):
    """ Error thrown when loading fails"""

def load_cif(filename: str | pathlib.Path,
             supercell: Supercell = TiledSupercell(),
             entry_index: int=0,
             verbose=False) -> Structure:
    """ Load a CIF file in as a structure """
    # Get the right bit of data

    cif = ReadCif(str(filename), grammar=None, permissive=True)

    names = [key for key in cif.keys()]

    if verbose:
        print("All entries:", names)

    for i in range(entry_index, len(names)):

        try:
            name = names[i]
        except IndexError:
            raise ValueError(f"Entry index out of range (file has {len(names)} entries)")

        if verbose:
            print(f"Trying entry {name}")

        data = cif[name]

        # Get the spacegroup and lattice parameters
        for spacegroup_key in ['_space_group_name_h-m_alt',
                                '_symmetry_space_group_name_h-m',
                                '_space_group.Patterson_name_h-m',
                                '_space_group.patterson_name_h-m']:

            try:
                spacegroup_name = data[spacegroup_key]
                break
            except KeyError:
                if verbose:
                    print(f"No entry: {spacegroup_key}")
        else:
            if verbose:
                print(f"Failed to find space group entry in '{name}'")
            continue

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

        if "_atom_type_radius_bond" in data and "_atom_type_symbol" in data:
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

    raise LoadError("No entries in CIF file")

def load_mcif(filename: str | pathlib.Path, entry_index: int=0, verbose: bool=False) -> Structure:

    """ Load a mCIF file in as a structure """

    cif = ReadCif(str(filename), grammar=None, permissive=True)

    names = [key for key in cif.keys()]

    if verbose:
        print("All entries:", names)

    for i in range(entry_index, len(names)):

        try:
            name = names[i]
        except IndexError:
            raise ValueError(f"Entry index out of range (file has {len(names)} entries)")

        if verbose:
            print(f"Trying entry {name}")

        data = cif[name]

        # Get the magnetic spacegroup

        operation_strings = data["_space_group_symop_magn_operation.xyz"]
        centering_strings = data["_space_group_symop_magn_centering.xyz"]

        if verbose:
            print("Magnetic Operations:")
            print(operation_strings)
            print("Centerings:")
            print(centering_strings)

        base_operations = [MagneticOperation.from_text(text) for text in operation_strings]
        centerings = [MagneticOperation.from_text(text) for text in centering_strings]

        sg = MagneticSpaceGroup.closure(base_operations, centerings)

        # Get the unit cell

        a = parse_float_entry(data["_cell_length_a"])
        b = parse_float_entry(data["_cell_length_b"])
        c = parse_float_entry(data["_cell_length_c"])

        alpha = parse_float_entry(data["_cell_angle_alpha"])
        beta = parse_float_entry(data["_cell_angle_beta"])
        gamma = parse_float_entry(data["_cell_angle_gamma"])

        cell = UnitCell(a,b,c,alpha=alpha,beta=beta,gamma=gamma)

        # Atom radii if given

        radius_lookup = None

        if "_atom_type_radius_bond" in data and "_atom_type_symbol" in data:
            radius_lookup = {atom: parse_float_entry(radius)
             for atom, radius in zip(data["_atom_type_symbol"], data["_atom_type_radius_bond"])}

        # Supercell components
        n_components = 1 # TODO: Actually read from file

        # Moments and positions

        sites = []
        moments = {}
        try:
            # Get the moments, using components in each crystal axis
            for label, mx, my, mz in zip(
                data["_atom_site_moment.label"],
                data["_atom_site_moment.crystalaxis_x"],
                data["_atom_site_moment.crystalaxis_y"],
                data["_atom_site_moment.crystalaxis_z"]):

                crystal_moment = np.array([
                    parse_float_entry(mx),
                    parse_float_entry(my),
                    parse_float_entry(mz)], dtype=float)

                # The following essentially orthogonalises the description
                ## Convert to unit cell
                # print("Lattice matrix:", cell._normalised_xyz)
                # print("check transpose:", cell._xyz)
                # print("Lattice matrix inv:", cell._normalised_xyz_inv)

                moment_ijk = cell._normalised_xyz @ crystal_moment

                # print("Moment in lattice:", moment_ijk)

                ## Convert spinw type coordinates
                moment = cell.spin_lattice_units_to_cartesian(moment_ijk)
                # print("Moment in xyz:", moment)
                #
                # print("Moment size:", np.sqrt(np.sum(moment**2)))

                moments[label] = moment

        except Exception as e:
            raise ValueError("Failed to read _atom_site_moment.label or"
                             "_atom_site_moment.crystalaxis_x/y/z, maybe"
                             "moment is specified in a different coordinate"
                             "system")

        # Get fourier components for multi-k structures

        # Create the lattice sites
        for label, atom, x_string, y_string, z_string in zip(
                data["_atom_site_label"],
                data["_atom_site_type_symbol"],
                data["_atom_site_fract_x"],
                data["_atom_site_fract_y"],
                data["_atom_site_fract_z"]):

            metadata = SiteMetadata.metadata_from_name(label)

            if radius_lookup is not None:
                metadata.radius = radius_lookup[atom]

            x = parse_float_entry(x_string) % 1
            y = parse_float_entry(y_string) % 1
            z = parse_float_entry(z_string) % 1

            if label in moments:
                moment = moments[label]
            else:
                moment = np.zeros((n_components, 3), dtype=float)

            sites.append(LatticeSite(x,y,z,
                                     supercell_spins=moment,
                                     name=label,
                                     metadata=metadata))



        return Structure(sites, unit_cell = cell, spacegroup=sg)

    raise LoadError("No entries in CIF file")

if __name__ == "__main__":
    structure = load_mcif("../examples/data/mcif/1.2_CuSe2O5.mcif", verbose=True)

    structure.print_summary()

    from pyspinw import view

    view(structure)