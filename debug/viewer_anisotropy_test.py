from multiprocessing import freeze_support

from pyspinw import AxisMagnitudeAnisotropy
from pyspinw.interface import generate_exchanges, axis_anisotropies
from pyspinw.exchange import HeisenbergExchange
from pyspinw.gui.viewer import show_object
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.site import LatticeSite
from pyspinw.structures import Structure
from pyspinw.symmetry.supercell import TiledSupercell
from pyspinw.symmetry.unitcell import UnitCell

if __name__ == "__main__":
    freeze_support()

    unit_cell = UnitCell(1,1,1)

    x = LatticeSite(0, 0.5, 0.5, 0, 0, 1, name="X")

    sites = [x]

    s = Structure(sites, unit_cell=unit_cell, supercell=TiledSupercell(scaling=(2, 2, 2)))

    anisotropies = [AxisMagnitudeAnisotropy(site=x, direction=[1,0,0], a=4)]

    exchanges = []
    hamiltonian = Hamiltonian(s, exchanges, anisotropies=anisotropies)

    print(hamiltonian.anisotropies)

    hamiltonian.print_summary()


    show_object(hamiltonian)
