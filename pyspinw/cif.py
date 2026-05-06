""" CIF File Loading """

from CifFile import ReadCif

from pyspinw.symmetry.group import NoSuchGroup
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.interface import spacegroup


def parse_float_entry(s: str):
    """ Parse a float entry, they can have brackets for values that go beyond the specified precision"""

    s = s.replace("(", "")
    s = s.replace(")", "")

    return float(s)

def load_cif(filename: str, entry_index=0):
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

    print(cell)
    print(sg)


if __name__ == "__main__":
    load_cif("../example_structures/1100231.cif")

