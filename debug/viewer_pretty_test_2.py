from multiprocessing import freeze_support

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

    unit_cell = UnitCell(1,1,1)

    x = LatticeSite(0, 0, 0, 0, 0, 1, name="X")
    y = LatticeSite(1, 0, 0, 0, 0, 1, name="Y")
    z = LatticeSite(0, 1, 0, 0, 0, 1, name="Z")
    w = LatticeSite(0,0,1,0,0,1, name="W")

    sites = [x, y, z, w]

    s = Structure(sites, unit_cell=unit_cell, supercell=TrivialSupercell(scaling=(1,1,1)))

    exchanges = couplings(sites=sites,
                           unit_cell=unit_cell,
                           max_distance=1.0001,
                           coupling_type=HeisenbergCoupling,
                           j=-1)

    anisotropies = axis_anisotropies([x], 1, (0,0,1)) + axis_anisotropies([x,y], 1, (0,1,0))

    hamiltonian = Hamiltonian(s, exchanges, anisotropies=anisotropies)

    hamiltonian.print_summary()


    show_hamiltonian(hamiltonian)