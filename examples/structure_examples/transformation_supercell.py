from pyspinw import *

supercell = TransformationSupercell([
        (CommensuratePropagationVector(1/3, 0, 0), RotationTransform([0, 0, 1])),
        (CommensuratePropagationVector(0, 1/2, 0), RotationTransform([0, 0, 1]))])


structure = Structure(
    [LatticeSite(0.5, 0.5, 0.5, 1, 0, 0)],
    UnitCell(1, 1, 1),
    spacegroup("p1"),
    supercell=supercell)

hamiltonian = Hamiltonian(structure, [])

view(hamiltonian)