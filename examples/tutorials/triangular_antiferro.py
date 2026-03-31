""" Triangular Antiferromagnet example """

from multiprocessing.spawn import freeze_support

from pyspinw.coupling import HeisenbergCoupling
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.interface import generate_exchanges, axis_anisotropies, generate_helical_structure
from pyspinw.path import Path
from pyspinw.symmetry.unitcell import UnitCell

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

    unit_cell = UnitCell(3, 3, 4, gamma=120)

    sites = generate_helical_structure(unit_cell, positions=[[0,0,0]], moments=[[0,1,0]], magnitudes=[3./2], names=['X'],
                                       perpendicular=[0,0,1], propagation_vector=[1./3., 1./3., 0])

    exchanges = generate_exchanges(sites=sites,
                                   bond=1,
                                   coupling_type=HeisenbergCoupling,
                                   j=1)

    anisotropies = axis_anisotropies(sites, 0.2)
    hamiltonian = Hamiltonian(sites, exchanges, anisotropies)

    hamiltonian.print_summary()

    path = Path([[0,0,0], [1,1,0]], n_points_per_segment=401)
    import matplotlib.pyplot as plt
    fig = hamiltonian.spaghetti_plot(path, show=False)
    fig.axes[0].set_ylim(0, 4)
    fig.axes[1].set_ylim(0, 3)
    plt.show()
