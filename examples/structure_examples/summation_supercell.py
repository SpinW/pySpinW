import numpy as np
from pyspinw import *

x = LatticeSite(0.5, 0.5, 0.5, supercell_spins=[[1,0,0], [0,1,0]])

supercell = SummationSupercell([
        CommensuratePropagationVector(1/3, 0, 0),
        CommensuratePropagationVector(1/3, 0, 0, phase=np.pi/2)])


structure = Structure([x], UnitCell(1, 1, 1), supercell=supercell)

hamiltonian = Hamiltonian(structure, [])

view(hamiltonian)