""" Kagome Ferromagnet example """

from multiprocessing.spawn import freeze_support

from pyspinw.coupling import HeisenbergCoupling
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

    unit_cell = UnitCell(1,1,1, gamma=60)
    group = spacegroup("p1")

    x = LatticeSite(0, 0, 0, 0, 0, 1, name="X")
    y = LatticeSite(0.5, 0, 0, 0, 0, 1, name="Y")
    z = LatticeSite(0, 0.5, 0, 0, 0, 1, name="Z")

    sites = [x, y, z]

    s = Structure(sites,
        unit_cell=unit_cell,
        spacegroup=group,
        supercell=TrivialSupercell(scaling=(3,3,1)))

    exchanges = []
    exchanges += couplings(sites=[x,y],
                          unit_cell=unit_cell,
                          max_distance=0.6,
                          coupling_type=HeisenbergCoupling,
                          j=-1)

    exchanges += couplings(sites=[y, z],
                           unit_cell=unit_cell,
                           max_distance=0.6,
                           coupling_type=HeisenbergCoupling,
                           j=-1)

    exchanges += couplings(sites=[z, x],
                           unit_cell=unit_cell,
                           max_distance=0.6,
                           coupling_type=HeisenbergCoupling,
                           j=-1)

    exchanges += couplings(sites=[y, x],
                           unit_cell=unit_cell,
                           max_distance=0.6,
                           coupling_type=HeisenbergCoupling,
                           j=-1)

    exchanges += couplings(sites=[z, y],
                           unit_cell=unit_cell,
                           max_distance=0.6,
                           coupling_type=HeisenbergCoupling,
                           j=-1)

    exchanges += couplings(sites=[x, z],
                           unit_cell=unit_cell,
                           max_distance=0.6,
                           coupling_type=HeisenbergCoupling,
                           j=-1)

    debug_plot(s, exchanges, show=False)

    hamiltonian = Hamiltonian(s, exchanges)

    path = Path([[-0.5,0,0], [0,0,0], [0,0.5,0.5]])

    print("Sites:")
    for site in sites:
        print(site)

    print("Couplings:")
    for exchange in exchanges:
        print(exchange, " vector =",exchange.vector(unit_cell=unit_cell))


    hamiltonian.energy_plot(path)
