""" Antiferromagnetic chain example """

from multiprocessing.spawn import freeze_support

from pyspinw.anisotropy import AxisMagnitudeAnisotropy
from pyspinw.coupling import HeisenbergCoupling
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.interface import spacegroup, couplings, filter, axis_anisotropies, axis_anisotropies
from pyspinw.path import Path
from pyspinw.site import LatticeSite
from pyspinw.symmetry.supercell import TrivialSupercell
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structures import Structure

from pyspinw.debug_plot import debug_plot

if __name__ == "__main__":
    freeze_support()

    unit_cell = UnitCell(2,1,1)

    sites = [LatticeSite(0, 0, 0, 0,0,1, name="X"),
             LatticeSite(0.5,0,0, 0,0, -1, name="Y")]

    s = Structure(sites, unit_cell=unit_cell)


    exchanges = couplings(sites=sites,
                          unit_cell=unit_cell,
                          max_distance=1.1,
                          coupling_type=HeisenbergCoupling,
                          j=1,
                          direction_filter=filter([1,0,0], symmetric=True))

    anisotropies = axis_anisotropies(sites, -0.1)

    hamiltonian = Hamiltonian(s, exchanges, anisotropies)

    hamiltonian.print_summary()

    path = Path([[0,0,0], [1,0,0]])
    hamiltonian.energy_plot(path, field=[0,0,7])
