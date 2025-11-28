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

    unit_cell = UnitCell(2,1,1)
    group = spacegroup("p1")

    sites = [LatticeSite(0, 0, 0, 0,0,1, name="X"),
             LatticeSite(0.5,0,0, 0,0, -1, name="Y")]

    s = Structure(sites,
        unit_cell=unit_cell,
        spacegroup=group,
        supercell=TrivialSupercell())


    exchanges = []
    exchanges += couplings(sites=sites,
                           unit_cell=unit_cell,
                           max_distance=1.1,
                           coupling_type=HeisenbergCoupling,
                           j=1,
                           direction_filter=filter([1,0,0], symmetric=True))

    exchanges += couplings(sites=sites[::-1],
                           unit_cell=unit_cell,
                           max_distance=1.1,
                           coupling_type=HeisenbergCoupling,
                           j=1,
                           direction_filter=filter([1, 0, 0], symmetric=True))

    print("Sites:")
    for site in sites:
        print(site)

    print("Couplings:")
    for exchange in exchanges:
        print(exchange, " vector =",exchange.vector(unit_cell=unit_cell))

    hamiltonian = Hamiltonian(s, exchanges)

    path = Path([[0,0,0], [1,0,0]])

    hamiltonian.energy_plot(path)