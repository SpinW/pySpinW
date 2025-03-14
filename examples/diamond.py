from ase.lattice import FCC

from ase.atoms import Atom

from pyspinw.sample import SingleCrystal
from pyspinw.structure import MagneticLattice
from pyspinw.hamiltonian import HeisenbergMagnet
from pyspinw.experiment import Experiment


lattice = FCC()
atoms = [Atom("C", (0,0,0)),
         Atom("C", (1/4, 1/4, 1/4))]

magnetic_structure = MagneticLattice(lattice, atoms, (2, 2, 2))
hamiltonian = HeisenbergMagnet(lattice=lattice, magnetic_structure=magnetic_structure, coupling=[])
sample = SingleCrystal(hamiltonian=hamiltonian)

experiment = Experiment(sample)

