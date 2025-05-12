from dataclasses import dataclass

from pyspinw.symmetry.system import CrystalSystem


class LatticeType:
    name: str
    letter: str

@dataclass
class Bravais:
    system: CrystalSystem
    lattice_type: LatticeType

lattice_types: list[LatticeType] = []