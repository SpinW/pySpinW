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

    x = LatticeSite(0.5, 0, 0, 1, 0, 0, name="X")
    y = LatticeSite(0, 0.5, 0, 0, 1, 0, name="Y")
    z = LatticeSite(0, 0, 0.5, 0, 0, 1, name="Z")

    sites = [x, y, z]

    s = Structure(sites, unit_cell=unit_cell)

    hamiltonian = Hamiltonian(s, [])

    hamiltonian.print_summary()


    show_hamiltonian(hamiltonian)