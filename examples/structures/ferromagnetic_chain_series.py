""" Ferromagnetic chain example """

from multiprocessing.spawn import freeze_support

from pyspinw.coupling import HeisenbergCoupling
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.interface import generate_exchanges, filter
from pyspinw.path import Path
from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structures import Structure
import sys

if __name__ == "__main__":
    freeze_support()

    use_rust = "py" not in sys.argv[1] if len(sys.argv) > 1 else True

    unit_cell = UnitCell(1,1,1)

    only_site = LatticeSite(0, 0, 0, 0,0,1, name="X")

    s = Structure([only_site], unit_cell=unit_cell)

    exchanges = generate_exchanges(sites=[only_site],
                                   unit_cell=unit_cell,
                                   max_distance=1.1,
                                   coupling_type=HeisenbergCoupling,
                                   j=-1,
                                   direction_filter=filter([1,0,0]))

    hamiltonian = Hamiltonian(s, exchanges)

    path = Path([[0,0,0], [1,0,0]])

    parameterized_hamiltonian = hamiltonian.parameterize((exchanges[0], "j"))

    for j in [-0.5, -1.0, -1.5]:
        ham = parameterized_hamiltonian(j)

        ham.print_summary()

    parameterized_hamiltonian.energy_plot([-0.5, -1.0, -1.5], path)

