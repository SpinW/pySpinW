""" Kagome Antiferromagnet example """

from multiprocessing.spawn import freeze_support

from pyspinw.coupling import HeisenbergCoupling
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.interface import spacegroup, couplings, filter
from pyspinw.path import Path
from pyspinw.site import LatticeSite
from pyspinw.symmetry.supercell import SummationSupercell, CommensuratePropagationVector, TrivialSupercell
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structures import Structure
from numpy import sqrt

from pyspinw.debug_plot import debug_plot

"""
AFkagome = spinw;
AFkagome.genlattice('lat_const',[6 6 10],'angled',[90 90 120],'spgr','P -3')
AFkagome.addatom('r',[1/2 0 0],'S', 1,'label','MCu1','color','r')
AFkagome.gencoupling('maxDistance',7)
AFkagome.addmatrix('label','J1','value',1.00,'color','r')
AFkagome.addmatrix('label','J2','value',0.11,'color','g')
AFkagome.addcoupling('mat','J1','bond',1)
AFkagome.addcoupling('mat','J2','bond',2)
S0 = [1 -2 1; 2 -1 -1; 0 0 0];
AFkagome.genmagstr('mode','direct','k',[0 0 0],'n',[0 0 1],'unit','lu','S',S0);
afkSpec = AFkagome.spinwave({[-1/2 0 0] [0 0 0] [1/2 1/2 0] 100},'hermit',true);
figure; subplot(211)
sw_plotspec(afkSpec,'mode',1,'colorbar',false,'colormap',[0 0 0],'dashed',true,'axLim',[0 3])
afkSpec = sw_egrid(sw_neutron(afkSpec),'Evect',linspace(0,6.5,500),'component','Sperp');
afkSpec = sw_omegasum(afkSpec,'zeroint',1e-6);
subplot(212); sw_plotspec(afkSpec,'mode',2,'log',false,'axLim',[0 3])
"""

if __name__ == "__main__":
    """Reproduces Tutorial 7: https://spinw.org/tutorials/07tutorial"""
    freeze_support()

    unit_cell = UnitCell(6, 6, 10, gamma=120)

    x = LatticeSite(0.5, 0, 0, 1, 2, 0, name="X", unit="lu")
    y = LatticeSite(0, 0.5, 0, -2, -1, 0, name="Y", unit="lu")
    z = LatticeSite(0.5, 0.5, 0, 1, -1, 0, name="Z", unit="lu")

    sites = [x, y, z]
    s = Structure(sites, unit_cell=unit_cell, supercell=TrivialSupercell(scaling=(1,1,1)))

    j1 = couplings(sites=[x, y, z], unit_cell=unit_cell, min_distance=0, max_distance=4.1, coupling_type=HeisenbergCoupling, j=1)
    j2 = couplings(sites=[x, y, z], unit_cell=unit_cell, min_distance=5, max_distance=5.2, coupling_type=HeisenbergCoupling, j=0.11)
    exchanges = j1 + j2

    debug_plot(s, exchanges, show=False)

    hamiltonian = Hamiltonian(s, exchanges)

    hamiltonian.print_summary()
    hamiltonian.expand()

    path = Path([[-0.5,0,0], [0,0,0], [0.5,0.5,0]])
    hamiltonian.plot(path)
