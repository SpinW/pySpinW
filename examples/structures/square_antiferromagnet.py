""" Square-lattice Antiferromagnet example """

from multiprocessing.spawn import freeze_support

from pyspinw.exchange import HeisenbergExchange
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.interface import generate_exchanges, generate_structure
from pyspinw.path import Path
from pyspinw.site import LatticeSite
from pyspinw.symmetry.supercell import SummationSupercell, CommensuratePropagationVector
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structures import Structure
import sys

"""
[J, Jp, Jpp, Jc] = deal(138.3, 2, 2, 38);
lacuo = sw_model('squareAF',[J-Jc/2 Jp-Jc/4 Jpp]/2,0);
lacuo.unit_cell.S = 1/2;
lacuoSpec = lacuo.spinwave({[3/4 1/4 0] [1/2 1/2 0] [1/2 0 0] [3/4 1/4 0] [1 0 0] [1/2 0 0] 100},'hermit',false);
lacuoSpec = sw_egrid(sw_neutron(lacuoSpec),'component','Sperp');
figure; subplot(2,1,1)
sw_plotspec(lacuoSpec,'mode','disp','axLim',[0 350],'dE',35,'dashed',true)
subplot(2,1,2)
lacuoSpec = sw_omegasum(lacuoSpec,'zeroint',1e-5,'tol',1e-3);
sw_plotspec(lacuoSpec,'mode',2,'axLim',[0 20],'dashed',true,'colormap',[0 0 0])
"""

if __name__ == "__main__":
    """Reproduces Tutorial 11: https://spinw.org/tutorials/11tutorial"""
    use_rust = "py" not in sys.argv[1] if len(sys.argv) > 1 else True

    unit_cell = UnitCell(3, 3, 9)

    sites = generate_structure(unit_cell, positions=[[0,0,0]], spins=[[1, 0, 0]], names=['X'],
                               propagation_vectors=[[0.5, 0.5, 0]])

    j1 = generate_exchanges(sites, bond=1, j=59.65)
    j2 = generate_exchanges(sites, bond=2, j=-3.75)
    j3 = generate_exchanges(sites, bond=3, j=1)
    exchanges = j1 + j2 + j3

    hamiltonian = Hamiltonian(sites, exchanges)
    hamiltonian.print_summary()

    path = Path([[3/4, 1/4, 0], [1/2, 1/2, 0], [1/2, 0, 0], [3/4, 1/4, 0], [1, 0, 0], [1/2, 0, 0]], n_points_per_segment=51)
    import matplotlib.pyplot as plt
    fig = hamiltonian.spaghetti_plot(path, show=False, use_rust=use_rust)
    fig.axes[0].set_ylim(0, 300)
    fig.axes[1].set_ylim(0, 0.5)
    plt.show()
