from multiprocessing.spawn import freeze_support

from pyspinw.hamiltonian import Hamiltonian
from pyspinw.coupling import HeisenbergCoupling
from pyspinw.couplinggroup import CouplingGroup
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.interface import spacegroup, couplings, filter
from pyspinw.path import Path
from pyspinw.site import LatticeSite
from pyspinw.symmetry.supercell import TrivialSupercell
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structures import Structure

from pyspinw.debug_plot import debug_plot

if __name__ == "__main__":
    freeze_support()

    unit_cell = UnitCell(1,1,1)
    group = spacegroup("p1")

    only_site = LatticeSite(0, 0, 0, 0,0,1, name="X")

    s = Structure([only_site],
        unit_cell=unit_cell,
        spacegroup=group,
        supercell=TrivialSupercell(scaling=(3,3,3)))


    exchanges = couplings(sites=[only_site],
                          unit_cell=unit_cell,
                          max_distance=1.1,
                          coupling_type=HeisenbergCoupling,
                          j=-1,
                          direction_filter=filter([1,0,0]))

    hamiltonian = Hamiltonian(s, exchanges)

    path = Path([[0,0,0], [1,0,0]])

    hamiltonian.energy_plot(path)