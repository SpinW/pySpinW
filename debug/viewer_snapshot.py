from multiprocessing import freeze_support

from pyspinw.interface import generate_exchanges, axis_anisotropies
from pyspinw.exchange import HeisenbergExchange
from pyspinw.gui.viewer import snapshot
from pyspinw.gui.displayoptions import DisplayOptions
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.site import LatticeSite
from pyspinw.structure import Structure
from pyspinw.symmetry.supercell import TiledSupercell
from pyspinw.symmetry.unitcell import UnitCell

if __name__ == "__main__":
    freeze_support()

    unit_cell = UnitCell(1,1,1, gamma=60)

    x = LatticeSite(0, 0, 0, 0, 0, 1, name="X")
    y = LatticeSite(0.5, 0, 0, 0, 0, 1, name="Y")
    z = LatticeSite(0, 0.5, 0, 0, 0, 1, name="Z")

    sites = [x, y, z]

    s = Structure(sites, unit_cell=unit_cell, supercell=TiledSupercell(scaling=(3, 3, 1)))

    exchanges = generate_exchanges(sites=[x, y, z],
                                   unit_cell=unit_cell,
                                   max_distance=0.6,
                                   exchange_type=HeisenbergExchange,
                                   j=-1)

    anisotropies = axis_anisotropies([x], 1, (0,0,1)) + axis_anisotropies([x,y], 1, (0,1,0))

    hamiltonian = Hamiltonian(s, exchanges, anisotropies=anisotropies)

    hamiltonian.print_summary()

    snapshot(hamiltonian,
             filename="snapshot_test.png",
             view_point=(0, 0, 10),
             display_options=DisplayOptions(show_sites=False))

    snapshot(hamiltonian, view_point=(10,0,0), display_options=DisplayOptions(default_exchange_color=(1,0,0), show_sites=False))

    import matplotlib.pyplot as plt

    plt.show()