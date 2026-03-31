from multiprocessing import freeze_support

import numpy as np

from pyspinw.calculations.energy_minimisation import ClassicalEnergyMinimisation, Free, Fixed
from pyspinw.interface import generate_exchanges
from pyspinw.coupling import HeisenbergCoupling
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.site import LatticeSite
from pyspinw.structures import Structure
from pyspinw.symmetry.supercell import TrivialSupercell
from pyspinw.symmetry.unitcell import UnitCell

if __name__ == "__main__":
    freeze_support()

    unit_cell = UnitCell(1,1,1, gamma=60)

    x1 = LatticeSite(0, 0, 0.5, 0, 0, 1, name="X1")
    x2 = LatticeSite(0, 0, 0, 0, 1, 0, name="X2")

    sites = [x1, x2]

    s = Structure(sites, unit_cell=unit_cell, supercell=TrivialSupercell())

    exchanges = generate_exchanges(sites=sites,
                                   unit_cell=unit_cell,
                                   max_distance=0.6,
                                   coupling_type=HeisenbergCoupling,
                                   j=1)

    hamiltonian = Hamiltonian(s, exchanges)

    hamiltonian.print_summary()

    minimiser = ClassicalEnergyMinimisation(hamiltonian, constraints=[Fixed, Free], field=np.array([0.0,0.0,0.0]))

    moment_history = [minimiser.moment_data.copy()]

    for i in range(100):
        # print(minimiser.moments)
        minimiser.iterate(0.1)
        print(minimiser.energy())

        moment_history.append(minimiser.moment_data.copy())

    print(minimiser.moment_data)

    import matplotlib.pyplot as plt
    ax = plt.figure().add_subplot(projection='3d')

    for i in range(len(sites)):
        xyz = np.squeeze(np.array([moments[i,:] for moments in moment_history]))
        ax.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], "--x")

    plt.show()
