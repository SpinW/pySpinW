from ase.lattice import FCC

from ase.atoms import Atom

from pyspinw.sample import SingleCrystal

lattice = FCC()
atoms = []

sample = SingleCrystal(atoms, lattice)