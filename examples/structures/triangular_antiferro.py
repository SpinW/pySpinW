""" Triangular Antiferromagnet example """

from multiprocessing.spawn import freeze_support

from pyspinw.coupling import HeisenbergCoupling
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.interface import spacegroup, couplings, filter
from pyspinw.interface import spacegroup, couplings, filter, axis_anisotropies
from pyspinw.path import Path
from pyspinw.site import LatticeSite
from pyspinw.symmetry.supercell import SummationSupercell, CommensuratePropagationVector
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structures import Structure
from math import sqrt

from pyspinw.debug_plot import debug_plot

"""
tri = spinw;
tri.genlattice('lat_const',[3 3 4],'angled',[90 90 120])
tri.addatom('r',[0 0 0],'S',3/2,'label','MCr3','color','orange')
tri.gencoupling()
tri.addmatrix('value',1,'label','J','color','SteelBlue')
tri.addcoupling('mat','J','bond',1)
tri.addmatrix('value',diag([0 0 0.2]),'label','D','color','r')
tri.addaniso('D')
tri.genmagstr('mode','helical','S',[0; 1; 0],'k',[1/3 1/3 0],'n', [0 0 1],'nExt',0.1);
triSpec = sw_neutron(tri.spinwave({[0 0 0] [1 1 0] 500}));
figure; subplot(211)
sw_plotspec(triSpec,'mode','disp','axLim',[0 7],'colormap',[0 0 0],'colorbar',false)
triSpec = sw_egrid(triSpec,'Evect',linspace(0,6.5,500),'component','Sperp');
triSpec = sw_omegasum(triSpec,'zeroint',1e-6);
subplot(212); sw_plotspec(triSpec,'mode',2,'log',false,'axLim',[0 3])
"""

if __name__ == "__main__":
    """Reproduces Tutorial 12: https://spinw.org/tutorials/12tutorial"""
    freeze_support()

    unit_cell = UnitCell(3, 3, 4, gamma=120)

    sites = [LatticeSite(0, 0, 0, -1j, 1, 0, name="X")]
    k = CommensuratePropagationVector(1./3., 1./3., 0)
    s = Structure(sites, unit_cell=unit_cell, supercell=SummationSupercell(propagation_vectors=[k]))

    exchanges = couplings(sites=sites,
                          unit_cell=unit_cell,
                          max_distance=3.1,
                          coupling_type=HeisenbergCoupling,
                          j=1)

    debug_plot(s, exchanges, show=False)

    anisotropies = axis_anisotropies(sites, 0.2)
    hamiltonian = Hamiltonian(s, exchanges, anisotropies)

    hamiltonian.print_summary()

    path = Path([[0,0,0], [1,1,0]], n_points_per_segment=401)
    import matplotlib.pyplot as plt
    fig = hamiltonian.spaghetti_plot(path, show=False)
    fig.axes[1].set_ylim(0, 5)
    plt.show()
