from dataclasses import dataclass

@dataclass
class LatticeType:
    name: str
    letter: str

PRIMITIVE = LatticeType(name="Primitive", letter="P")
BASE_CENTERED = LatticeType(name="Base Centered", letter="S")
BODY_CENTERED = LatticeType(name="Body Centered", letter="I")
FACE_CENTERED = LatticeType(name="Face Centered", letter="F")
RHOMBOHEDRAL = LatticeType(name="Rhombohedral", letter="R")

lattice_types: list[LatticeType] = [
    PRIMITIVE,
    BASE_CENTERED,
    BODY_CENTERED,
    FACE_CENTERED,
    RHOMBOHEDRAL
]

lattice_type_name_lookup = {bravais.name: bravais for bravais in lattice_types}
lattice_type_letter_lookup = {bravais.name: bravais for bravais in lattice_types}
