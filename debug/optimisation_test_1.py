from multiprocessing import freeze_support

import numpy as np

from pyspinw.calculations.optimisation.energy_minimisation import ClassicalEnergyMinimisation, Free
from pyspinw.interface import couplings, axis_anisotropies
from pyspinw.coupling import HeisenbergCoupling
from pyspinw.gui.viewer import show_hamiltonian
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.site import LatticeSite
from pyspinw.structures import Structure
from pyspinw.symmetry.supercell import TrivialSupercell
from pyspinw.symmetry.unitcell import UnitCell

if __name__ == "__main__":
    freeze_support()

    unit_cell = UnitCell(1,1,1, gamma=60)

    x1 = LatticeSite(0, 0, 0.5, 0, 0, 1, name="X1")
    x2 = LatticeSite(0, 0, 0, 0, 0, 1, name="X2")

    sites = [x1, x2]

    s = Structure(sites, unit_cell=unit_cell, supercell=TrivialSupercell())

    exchanges = couplings(sites=sites,
                           unit_cell=unit_cell,
                           max_distance=0.6,
                           coupling_type=HeisenbergCoupling,
                           j=1)

    hamiltonian = Hamiltonian(s, exchanges)

    hamiltonian.print_summary()

    minimiser = ClassicalEnergyMinimisation(hamiltonian, constraints=[Free, Free], field=np.array([0.0,0.0,0.0]))

    minimiser.iterate(0.01)