from multiprocessing import freeze_support

from pyspinw import AxisMagnitudeAnisotropy, Anisotropy
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

    x = LatticeSite(0.5, 0.75, 0.75, 0, 0, 1, name="S1")
    y = LatticeSite(0.5, 0.25, 0.75, 0, 0, 1, name="S2")
    z = LatticeSite(0.5, 0.75, 0.25, 0, 0, 1, name="S3")
    w = LatticeSite(0.5, 0.25, 0.25, 0, 0, 1, name="S4")

    sites = [x,y,z,w]

    s = Structure(sites, unit_cell=unit_cell, supercell=TiledSupercell(scaling=(1, 1, 1)))

    anisotropies = [
        AxisMagnitudeAnisotropy(site=x, direction=[1,0,0], a=1),
        AxisMagnitudeAnisotropy(site=y, direction=[1,0,0], a=-1),
        Anisotropy(site=z, anisotropy_matrix=[[1,0,0],[-1,0,0],[0,0,0]]),
        AxisMagnitudeAnisotropy(site=w, direction=[1,0,0], a=0),
    ]

    exchanges = []
    hamiltonian = Hamiltonian(s, exchanges, anisotropies=anisotropies)

    print(hamiltonian.anisotropies)

    hamiltonian.print_summary()


    show_object(hamiltonian)
