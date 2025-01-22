from importlib import resources

from dataclasses import dataclass

from Tools.scripts.stable_abi import generators

""" Deals with magnetic symmetry group names and numbers """

@dataclass
class MagneticSymmetryGroup:
    """ Class containing the many different ways of notating the same symmetry group"""
    litvin: int
    bns_number: float
    bns_symbol: str
    og_number: str
    og_symbol: str
    int_number: int
    int_symbol: str
    mhall: str
    uni_symbol: str
    generators: str


class _MagneticSymmetryLookup:
    """ Lookup table for converting between types of notation """
    def __init__(self):
        self.litvin: dict[int, MagneticSymmetryGroup] = {}
        self.bns_number: dict[float, MagneticSymmetryGroup] = {}
        self.bns_symbol: dict[str, MagneticSymmetryGroup] = {}
        self.og_number: dict[str, MagneticSymmetryGroup] = {}
        self.og_symbol: dict[str, MagneticSymmetryGroup] = {}
        self.int_number: dict[int, MagneticSymmetryGroup] = {}
        self.int_symbol: dict[str, MagneticSymmetryGroup] = {}
        self.mhall: dict[str, MagneticSymmetryGroup] = {}
        self.uni_symbol: dict[str, MagneticSymmetryGroup] = {}

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
        self.uni_symbol[entry.uni_symbol] = entry


name_converter = _MagneticSymmetryLookup()

# Load the bns/uni data - taken from S3 of https://doi.org/10.1107/S2053273321012912
_bns_uni_lookup: dict[float, tuple[str, str, str]] = {}
with resources.open_text("pyspinw.data", "bns_uni_hall_gen.txt") as file:
    for line in file:
        try:
            bns = float(line[:13])
            uni = line[13:44].strip()
            hall = line[44:68].strip()
            gen = line[68:].strip()

            _bns_uni_lookup[bns] = (uni, hall, gen)

        except:
            pass

# Load the other list, make the objects, add to the name converter
with resources.open_text("pyspinw.data", "magnetic_symmetry_conventions.csv") as file:
    for line in file:
        parts = [part.strip() for part in line.split(",")]

        bns_number = float(parts[1])

        name_converter._add_entry(
            MagneticSymmetryGroup(
                litvin=int(parts[0]),
                bns_number=bns_number,
                bns_symbol=parts[2],
                og_number=parts[3],
                og_symbol=parts[4],
                int_number=int(parts[5]),
                int_symbol=parts[6],
                mhall=parts[7],
                uni_symbol=_bns_uni_lookup[bns_number][0],
                generators=_bns_uni_lookup[bns_number][2]
            )
        )


if __name__ == "__main__":
    # Example
    print(name_converter.bns_number[14.84].bns_symbol) # should be P_c 21/c
    print(name_converter.bns_number[14.84].litvin)