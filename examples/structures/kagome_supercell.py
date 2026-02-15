""" Kagome 3x3 Antiferromagnet example """

from multiprocessing.spawn import freeze_support

from pyspinw.coupling import HeisenbergCoupling
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.interface import spacegroup, couplings, filter
from pyspinw.path import Path
from pyspinw.site import LatticeSite
from pyspinw.symmetry.supercell import SummationSupercell, CommensuratePropagationVector
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structures import Structure
from math import sqrt

from pyspinw.debug_plot import debug_plot

"""
AF33kagome = spinw;
AF33kagome.genlattice('lat_const',[6 6 40],'angled',[90 90 120],'spgr','P -3')
AF33kagome.addatom('r',[1/2 0 0],'S', 1,'label','MCu1','color','r')
AF33kagome.gencoupling('maxDistance',7)
AF33kagome.addmatrix('label','J1','value',1.00,'color','g')
AF33kagome.addcoupling('mat','J1','bond',1)
S0 = [0 0 -1;
      1 1 -1;
      0 0 0];
AF33kagome.genmagstr('mode','helical','k',[-1/3 -1/3 0],...
    'n',[0 0 1],'unit','lu','S',S0,'nExt',0.1);
kag33Spec = AF33kagome.spinwave({[-1/2 0 0] [0 0 0] [1/2 1/2 0] 250},'hermit',false);
figure; sw_plotspec(kag33Spec)
"""

if __name__ == "__main__":
    """Reproduces Tutorial 8: https://spinw.org/tutorials/08tutorial"""
    freeze_support()

    unit_cell = UnitCell(6, 6, 40, gamma=120)

    x = LatticeSite(0.5, 0,   0, 0, 1, 0, name="X", unit="lu")
    y = LatticeSite(0,   0.5, 0, 0, 1, 0, name="Y", unit="lu")
    z = LatticeSite(0.5, 0.5, 0, -1, -1, 0, name="Z", unit="lu")

    sites = [x, y, z]
    k = CommensuratePropagationVector(-1./3., -1./3., 0)
    s = Structure(sites, unit_cell=unit_cell, supercell=SummationSupercell(propagation_vectors=[k]))

    exchanges = couplings(sites=[x, y, z],
                          unit_cell=unit_cell,
                          max_distance=3.1,
                          coupling_type=HeisenbergCoupling,
                          j=1)

    debug_plot(s, exchanges, show=False)

    hamiltonian = Hamiltonian(s, exchanges)

    hamiltonian.print_summary()

    path = Path([[-0.5,0,0], [0,0,0], [0.5,0.5,0]])
    hamiltonian.spaghetti_plot(path, use_rust=False)
