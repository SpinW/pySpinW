from typing import Callable
from dataclasses import dataclass

import numpy as np
from CifFile import ReadCif, CifFile
import spglib

@dataclass
class CellData:
    a: float
    b: float
    c: float
    alpha: float
    beta: float
    gamma: float


def extract_cell_data(data: CifFile) -> CellData:
    a = data["_cell_length_a"]
    b = data["_cell_length_b"]
    c = data["_cell_length_c"]
    alpha = data["_cell_angle_alpha"]
    beta = data["_cell_angle_beta"]
    gamma = data["_cell_angle_gamma"]

    return CellData(a, b, c, alpha, beta, gamma)

def load_mcif(filename: str, block_name: str | None = None):
    """ Create a unit cell based on the data in an mcif file

    :param filename: string with filename to be opened
    :param block_name: optional string, name of dataset within file, if not present, this method will give the first one
    """

    data: CifFile = ReadCif(filename)

    # Get the block
    if block_name is None:
        block = data.first_block()
    else:
        block = data[block_name]

    cell_data = extract_cell_data(block)

    print(cell_data)
    return data

def load_cif(filename: str):
    """ Create a unit cell based on the data in a cif file """


if __name__ == "__main__":
    data = load_mcif("../../example_structures/1.2_CuSe2O5.mcif")
    # data = load_mcif("../../example_structures/1.2_CuSe2O5-extended.mcif")
    # data = load_mcif("../../example_structures/1100231.cif")

    # print(data)