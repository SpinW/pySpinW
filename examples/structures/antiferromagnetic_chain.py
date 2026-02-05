""" Antiferromagnetic chain example """

from multiprocessing.spawn import freeze_support

from pyspinw.coupling import HeisenbergCoupling
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.interface import spacegroup, couplings, filter
from pyspinw.path import Path
from pyspinw.site import LatticeSite
from pyspinw.symmetry.supercell import SummationSupercell, CommensuratePropagationVector
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structures import Structure

from pyspinw.debug_plot import debug_plot

if __name__ == "__main__":
    """Reproduces Tutorial 2: https://spinw.org/tutorials/02tutorial"""
    freeze_support()

    unit_cell = UnitCell(3, 8, 8)

    sites = [LatticeSite(0, 0, 0, 0, 0, 1, name="MCu1")]

    k = CommensuratePropagationVector(0.5, 0, 0)
    s = Structure(sites, unit_cell=unit_cell, supercell=SummationSupercell(propagation_vectors=[k]))

    exchanges = couplings(sites=sites,
                          unit_cell=unit_cell,
                          max_distance=3.1,
                          coupling_type=HeisenbergCoupling,
                          j=1)


    hamiltonian = Hamiltonian(s, exchanges)
    hamiltonian.print_summary()

    path = Path([[0,0,0], [1,0,0]])
    hamiltonian.energy_plot(path, show_intensities=True)
