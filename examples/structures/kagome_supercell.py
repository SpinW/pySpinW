""" Kagome 3x3 Antiferromagnet example """

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
figure; subplot(211)
sw_plotspec(kag33Spec,'mode',1,'colorbar',false,'colormap',[0 0 0],'dashed',true,'axLim',[0 3])
kag33Spec = sw_egrid(sw_neutron(kag33Spec),'Evect',linspace(0,6.5,500),'component','Sperp');
kag33Spec = sw_omegasum(kag33Spec,'zeroint',1e-6);
subplot(212); sw_plotspec(kag33Spec,'mode',2,'log',false,'axLim',[0 3])
"""

if __name__ == "__main__":
    """Reproduces Tutorial 8: https://spinw.org/tutorials/08tutorial"""
    freeze_support()

    use_rust = "py" not in sys.argv[1] if len(sys.argv) > 1 else True

    unit_cell = UnitCell(6, 6, 40, gamma=120)

    x = LatticeSite(0.5, 0,   0, 0, 1, 0, S=1, name="X")
    y = LatticeSite(0,   0.5, 0, 0, 1, 0, S=1, name="Y")
    z = LatticeSite(0.5, 0.5, 0, -1, -1, 0, S=1, name="Z")
    s = genmagstr([x, y, z], unit_cell, mode='helical', k=[-1./3, -1./3, 0], n=[0, 0, 1], unit='lu')

    exchanges = couplings(sites=[x, y, z],
                          unit_cell=unit_cell,
                          max_distance=3.1,
                          coupling_type=HeisenbergCoupling,
                          j=1)

    hamiltonian = Hamiltonian(s, exchanges)

    hamiltonian.print_summary()

    #from pyspinw.gui.viewer import show_hamiltonian
    #show_hamiltonian(hamiltonian)

    path = Path([[-0.5,0,0], [0,0,0], [0.5,0.5,0]])
    import matplotlib.pyplot as plt
    fig = hamiltonian.spaghetti_plot(path, show=False, use_rust=use_rust)
    fig.axes[0].set_ylim(0, 3)
    fig.axes[1].set_ylim(0, 1)
    plt.show()
