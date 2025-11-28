from typing import Callable, Any
from dataclasses import dataclass

import numpy as np

import spglib

from CifFile import ReadCif, CifFile
from CifFile.StarFile import StarBlock



@dataclass
class CellData:
    a: float
    b: float
    c: float
    alpha: float
    beta: float
    gamma: float


@dataclass
class SymmetryData:
    space_group_number: int
    magnetic_group_number: int


@dataclass
class AtomTemp:
    name: str
    position: np.ndarray[Any, np.dtype[float]]
    magnetic_moment: np.ndarray[Any, np.dtype[float]] | None = None

def cif_float(string: str) -> float:
    """ Parse decimals in cif files - they might have brackets, because of course"""
    trimmed = string.replace("(","").replace(")","")
    return float(trimmed)

def extract_cell_data(block: StarBlock) -> CellData:
    """ Extract the data describing the unit cell from a cif file block"""
    a = block["_cell_length_a"]
    b = block["_cell_length_b"]
    c = block["_cell_length_c"]
    alpha = block["_cell_angle_alpha"]
    beta = block["_cell_angle_beta"]
    gamma = block["_cell_angle_gamma"]

    return CellData(a, b, c, alpha, beta, gamma)

def extract_site_data(block: StarBlock) -> list[AtomTemp]:
    """ Extract the data describing the atoms from a cif file block"""
    labels = block["_atom_site_label"]
    xs = block["_atom_site_fract_x"]
    ys = block["_atom_site_fract_y"]
    zs = block["_atom_site_fract_z"]

    atoms = {}
    for label, x, y, z in zip(labels, xs, ys, zs):

        position = np.array([
                    cif_float(x),
                    cif_float(y),
                    cif_float(z)])

        atoms[label] = AtomTemp(label, position)


    try:
        magnetic_labels = block["_atom_site_moment.label"]
        mxs = block["_atom_site_moment.crystalaxis_x"]
        mys = block["_atom_site_moment.crystalaxis_y"]
        mzs = block["_atom_site_moment.crystalaxis_z"]

        for label, mx, my, mz in zip(magnetic_labels, mxs, mys, mzs):
            m = np.array(
                [cif_float(mx),
                        cif_float(my),
                        cif_float(mz)], dtype=float)

            atoms[label].magnetic_moment = m

    except:
        pass

    return list(atoms.values())

def extract_symmetry(block: StarBlock) -> SymmetryData:
    pass

def load_mcif(filename: str, block_name: str | None = None):
    """ Create a unit cell based on the data in an mcif file

    :param filename: string with filename to be opened
    :param block_name: optional string, name of dataset within file, if not present, this method will give the first one
    """
    data: CifFile = ReadCif(filename)

    # Get the appropriate block
    if block_name is None:
        block = data.first_block()
    else:
        block = data[block_name]

    cell_data = extract_cell_data(block)
    atom_data = extract_site_data(block)

    # Apply the symmetry operations

    print(cell_data)
    print(atom_data)
    return data



def load_cif(filename: str):
    """ Create a unit cell based on the data in a cif file """


if __name__ == "__main__":
    data = load_mcif("../../example_structures/1.2_CuSe2O5.mcif")
    # data = load_mcif("../../example_structures/1.2_CuSe2O5-extended.mcif")
    # data = load_mcif("../../example_structures/1100231.cif")

    # print(data)
