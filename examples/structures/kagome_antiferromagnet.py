""" Kagome Antiferromagnet example """

from multiprocessing.spawn import freeze_support

from pyspinw.coupling import HeisenbergCoupling
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.interface import couplings
from pyspinw.path import Path
from pyspinw.site import LatticeSite
from pyspinw.legacy.genmagstr import genmagstr
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structures import Structure
import sys

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

    use_rust = "py" not in sys.argv[1] if len(sys.argv) > 1 else True

    unit_cell = UnitCell(6, 6, 10, gamma=120)

    x = LatticeSite(0.5, 0, 0, 1, 2, 0, S=1, name="X")
    y = LatticeSite(0, 0.5, 0, -2, -1, 0, S=1, name="Y")
    z = LatticeSite(0.5, 0.5, 0, 1, -1, 0, S=1, name="Z")
    s = genmagstr([x, y, z], unit_cell, mode='tile', unit='lu')

    j1 = couplings(sites=[x, y, z], unit_cell=unit_cell, min_distance=0, max_distance=4.1, coupling_type=HeisenbergCoupling, j=1)
    j2 = couplings(sites=[x, y, z], unit_cell=unit_cell, min_distance=5, max_distance=5.2, coupling_type=HeisenbergCoupling, j=0.11)
    exchanges = j1 + j2

    hamiltonian = Hamiltonian(s, exchanges)

    hamiltonian.print_summary()

    path = Path([[-0.5,0,0], [0,0,0], [0.5,0.5,0]])
    hamiltonian.spaghetti_plot(path, use_rust=use_rust)
