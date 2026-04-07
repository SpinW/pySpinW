from multiprocessing import freeze_support

from pyspinw.exchange import HeisenbergExchange
from pyspinw.gui.viewer import show_hamiltonian
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.site import LatticeSite
from pyspinw.structures import Structure
from pyspinw.symmetry.unitcell import UnitCell

if __name__ == "__main__":
    freeze_support()

    unit_cell = UnitCell(1,1,1, alpha=50, beta=60, gamma=70)

    origin = LatticeSite(0,0,0, name="origin")
    a = LatticeSite(1,0,0, name="a")
    b = LatticeSite(0,1,0, name="b")
    c = LatticeSite(0,0,1, name="c")

    sites = [origin, a, b, c]
    exchanges = [HeisenbergExchange(origin, a, j=0, name="a"),
                 HeisenbergExchange(origin, b, j=0, name="b"),
                 HeisenbergExchange(origin, c, j=0, name="c")]


    s = Structure(sites, unit_cell=unit_cell)

    hamiltonian = Hamiltonian(s, exchanges)

    hamiltonian.print_summary()


    show_hamiltonian(hamiltonian)