from multiprocessing import freeze_support

from pyspinw.interface import generate_exchanges, axis_anisotropies
from pyspinw.exchange import HeisenbergExchange
from pyspinw.gui.viewer import show_object
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.site import Site
from pyspinw.structure import Structure
from pyspinw.symmetry.supercell import TiledSupercell
from pyspinw.symmetry.unitcell import UnitCell

if __name__ == "__main__":
    freeze_support()

    unit_cell = UnitCell(1,1,1, gamma=60)

    x = Site(0, 0, 0, 0, 0, 1, name="X")
    y = Site(0.5, 0, 0, 0, 0, 1, name="Y")
    z = Site(0, 0.5, 0, 0, 0, 1, name="Z")
    w = Site(0.5, 0.5, 0, 0, 0, 0, name="W")

    sites = [x, y, z, w]

    s = Structure(sites, unit_cell=unit_cell, supercell=TiledSupercell(scaling=(3, 3, 1)))

    exchanges = generate_exchanges(sites=[x, y, z],
                                   unit_cell=unit_cell,
                                   max_distance=0.6,
                                   exchange_type=HeisenbergExchange,
                                   j=-1)

    hamiltonian = Hamiltonian(s, exchanges)

    hamiltonian.print_summary()


    show_object(hamiltonian)
