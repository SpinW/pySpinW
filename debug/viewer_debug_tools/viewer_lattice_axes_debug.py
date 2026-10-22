from multiprocessing import freeze_support

from pyspinw.exchange import HeisenbergExchange
from pyspinw.gui.viewer import show_object
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.site import Site
from pyspinw.structure import Structure
from pyspinw.symmetry.unitcell import UnitCell

if __name__ == "__main__":
    freeze_support()

    unit_cell = UnitCell(1,1,1, alpha=50, beta=60, gamma=70)

    origin = Site(0, 0, 0, name="origin")
    a = Site(1, 0, 0, name="a")
    b = Site(0, 1, 0, name="b")
    c = Site(0, 0, 1, name="c")

    sites = [origin, a, b, c]
    exchanges = [HeisenbergExchange(origin, a, j=0, name="a"),
                 HeisenbergExchange(origin, b, j=0, name="b"),
                 HeisenbergExchange(origin, c, j=0, name="c")]


    s = Structure(sites, unit_cell=unit_cell)

    hamiltonian = Hamiltonian(s, exchanges)

    hamiltonian.print_summary()


    show_object(hamiltonian)
