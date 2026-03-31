""" Triangular Antiferromagnet example """

from multiprocessing.spawn import freeze_support

from pyspinw.coupling import HeisenbergCoupling
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.interface import couplings, axis_anisotropies
from pyspinw.path import Path
from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structures import Structure
from pyspinw.legacy.genmagstr import genmagstr
import sys

"""
tri = spinw;
tri.genlattice('lat_const',[3 3 8],'angled',[90 90 120])
tri.addatom('r',[0 0 0],'S',3/2,'label','MCr3','color','orange')
tri.addatom('r',[0 0 0.5],'S',3/2,'label','MCr3','color','orange')
tri.gencoupling()
tri.addmatrix('value',1,'label','J','color','SteelBlue')
tri.addcoupling('mat','J','bond',1)
tri.addmatrix('value',-0.1,'label','J2','color','Cyan'); tri.addcoupling('mat','J2','bond',2)
tri.addmatrix('value',diag([0 0 0.2]),'label','D','color','r')
tri.addaniso('D')
tri.genmagstr('mode','helical','S',[0 0; 1 1; 0 0],'k',[1/3 1/3 0],'n', [0 0 1])
triSpec = sw_neutron(tri.spinwave({[0 0 0] [1 1 0] 400}));
figure; subplot(211)
sw_plotspec(triSpec,'mode','disp','axLim',[0 7],'colormap',[0 0 0],'colorbar',false)
triSpec = sw_egrid(triSpec,'Evect',linspace(0,6.5,500),'component','Sperp');
triSpec = sw_omegasum(triSpec,'zeroint',1e-6);
subplot(212); sw_plotspec(triSpec,'mode',2,'log',false,'axLim',[0 5])
"""

if __name__ == "__main__":
    """Reproduces Tutorial 12: https://spinw.org/tutorials/12tutorial"""

    unit_cell = UnitCell(3, 3, 8, gamma=120)

    sites = [LatticeSite(0, 0, 0, 0, 1, 0, name="X"), LatticeSite(0, 0, 0.5, 0, 1, 0, name="Y")]
    s = genmagstr(sites, unit_cell, magnitude=[3./2, 3./2],mode='helical', k=[1./3, 1./3, 0], n=[0, 0, 1])

    exchanges = couplings(sites=sites, unit_cell=unit_cell, max_distance=3.1, coupling_type=HeisenbergCoupling, j=1) \
              + couplings(sites=sites, unit_cell=unit_cell, min_distance=3.1, max_distance=4.1, j=-0.1, coupling_type=HeisenbergCoupling)

    anisotropies = axis_anisotropies(sites, 0.2)
    hamiltonian = Hamiltonian(s, exchanges, anisotropies)

    hamiltonian.print_summary()

    path = Path([[0,0,0], [1,1,0]], n_points_per_segment=401)
    import matplotlib.pyplot as plt
    fig = hamiltonian.spaghetti_plot(path, show=False, use_rotating=True)
    fig.axes[0].set_ylim(0, 7)
    fig.axes[1].set_ylim(0, 5)
    plt.show()
