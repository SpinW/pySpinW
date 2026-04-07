""" Triangular Antiferromagnet example """

from multiprocessing.spawn import freeze_support

from pyspinw.exchange import HeisenbergExchange
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.interface import generate_exchanges, generate_helical_structure, axis_anisotropies
from pyspinw.path import Path
from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structures import Structure
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
    use_rust = "py" not in sys.argv[1] if len(sys.argv) > 1 else True

    unit_cell = UnitCell(3, 3, 8, gamma=120)

    sites = generate_helical_structure(unit_cell, positions=[[0, 0, 0], [0, 0, 0.5]], moments=[[0, 1, 0], [0, 1, 0]],
                magnitudes=[3./2, 3./2], propagation_vector=[1./3, 1./3, 0], perpendicular=[0, 0, 1])

    exchanges = generate_exchanges(sites=sites, bond=1, exchange_type=HeisenbergExchange, j=1) \
              + generate_exchanges(sites=sites, bond=2, exchange_type=HeisenbergExchange, j=-0.1)

    anisotropies = axis_anisotropies(sites, 0.2)
    hamiltonian = Hamiltonian(sites, exchanges, anisotropies)

    hamiltonian.print_summary()

    path = Path([[0,0,0], [1,1,0]], n_points_per_segment=401)
    import matplotlib.pyplot as plt
    fig = hamiltonian.spaghetti_plot(path, show=False, use_rust=use_rust, use_rotating=True)
    fig.axes[0].set_ylim(0, 7)
    fig.axes[1].set_ylim(0, 5)
    plt.show()
