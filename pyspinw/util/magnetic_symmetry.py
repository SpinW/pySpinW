from importlib import resources

from dataclasses import dataclass

""" Deals with magnetic symmetry group names and numbers """

@dataclass
class MagneticSymmetryGroup:
    litvin: int
    bns_number: float
    bns_symbol: str
    og_number: str
    og_symbol: str
    int_number: int
    int_symbol: str
    mhall: str

class _MagneticSymmetryLookup:
    def __init__(self):
        self.litvin: dict[int, MagneticSymmetryGroup] = {}
        self.bns_number: dict[float, MagneticSymmetryGroup] = {}
        self.bns_symbol: dict[str, MagneticSymmetryGroup] = {}
        self.og_number: dict[str, MagneticSymmetryGroup] = {}
        self.og_symbol: dict[str, MagneticSymmetryGroup] = {}
        self.int_number: dict[int, MagneticSymmetryGroup] = {}
        self.int_symbol: dict[str, MagneticSymmetryGroup] = {}
        self.mhall: dict[str, MagneticSymmetryGroup] = {}

    def _add_entry(self, entry: MagneticSymmetryGroup):
        """ Add an entry to the lookup """
        self.litvin[entry.litvin] = entry
        self.bns_number[entry.bns_number] = entry
        self.bns_symbol[entry.bns_symbol] = entry
        self.og_number[entry.og_number] = entry
        self.og_symbol[entry.og_symbol] = entry
        self.int_number[entry.int_number] = entry
        self.int_symbol[entry.int_symbol] = entry
        self.mhall[entry.mhall] = entry


name_converter = _MagneticSymmetryLookup()

with resources.open_text("pyspinw.data", "magnetic_symmetry_conventions.csv") as file:
    for line in file:
        parts = [part.strip() for part in line.split(",")]

        name_converter._add_entry(
            MagneticSymmetryGroup(
                litvin=int(parts[0]),
                bns_number=float(parts[1]),
                bns_symbol=parts[2],
                og_number=parts[3],
                og_symbol=parts[4],
                int_number=int(parts[5]),
                int_symbol=parts[6],
                mhall=parts[7]
            )
        )


if __name__ == "__main__":
    # Example
    print(name_converter.bns_number[14.84].bns_symbol) # should be P_c 21/c
    print(name_converter.bns_number[14.84].litvin)