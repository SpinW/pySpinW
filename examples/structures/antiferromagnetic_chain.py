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

"""
AFMchain = spinw;
AFMchain.genlattice('lat_const',[3 8 8],'angled',[90 90 90],'spgr',0);
AFMchain.addatom('r',[0 0 0],'S',1,'label','MCu1','color','blue');
AFMchain.gencoupling('maxDistance',7)
AFMchain.addmatrix('label','Ja','value',1,'color','red');
AFMchain.addcoupling('mat','Ja','bond',1);
AFMchain.genmagstr('mode','direct','k',[1/2 0 0],'n',[1 0 0],'S',[0 0; 1 -1;0 0]);
afcSpec = AFMchain.spinwave({[0 0 0] [1 0 0] 523}, 'hermit',true);
figure; subplot(2,1,1)
sw_plotspec(afcSpec,'mode',4,'dE',0.2,'axLim',[0 3])
afcSpec = sw_egrid(sw_neutron(afcSpec),'Evect',linspace(0,6.5,500),'component','Sperp');
afcSpec = sw_omegasum(afcSpec,'zeroint',1e-6);
subplot(2,1,2)
sw_plotspec(afcSpec,'mode',2,'log',true,'axLim',[-4 10])
"""

if __name__ == "__main__":
    """Reproduces Tutorial 2: https://spinw.org/tutorials/02tutorial"""
    freeze_support()

    use_rust = "py" not in sys.argv[1] if len(sys.argv) > 1 else True

    unit_cell = UnitCell(3, 8, 8)

    sites = [LatticeSite(0, 0, 0, 0, 1, 0, name="MCu1")]

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
    hamiltonian.spaghetti_plot(path, scale='log', use_rust=use_rust)
