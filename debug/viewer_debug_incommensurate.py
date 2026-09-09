from multiprocessing import freeze_support

import numpy as np

from pyspinw import RotationSupercell
from pyspinw.exchange import HeisenbergExchange
from pyspinw.gui.viewer import show_object
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.site import LatticeSite
from pyspinw.structures import Structure
from pyspinw.symmetry.unitcell import UnitCell

unit_cell = UnitCell(1,1,1)

site = LatticeSite(0.5, 0.5, 0.5, 1, 0, 0, name="origin")

supercell = RotationSupercell(
    np.array([0,1,0]),
    np.array([0,0,np.sqrt(1/10)]))

s = Structure([site], unit_cell=unit_cell, supercell=supercell)

hamiltonian = Hamiltonian(s, [])

hamiltonian.print_summary()


show_object(hamiltonian)
