""" Antiferromagnetic chain example with applied magnetic field """

from multiprocessing.spawn import freeze_support

from pyspinw.anisotropy import AxisMagnitudeAnisotropy
from pyspinw.coupling import HeisenbergCoupling
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.interface import spacegroup, couplings, filter, axis_anisotropies, axis_anisotropies
from pyspinw.path import Path
from pyspinw.site import LatticeSite
from pyspinw.symmetry.supercell import TrivialSupercell
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structures import Structure

from pyspinw.debug_plot import debug_plot

"""
afc = spinw;
afc.genlattice('lat_const',[4 6 6])
afc.addatom('r',[  0 0 0],'S',1)
afc.addatom('r',[1/2 0 0],'S',1)
afc.addmatrix('label','J','value',1)
afc.gencoupling
afc.addcoupling('mat','J','bond',1)
afc.addmatrix('label','A','value',diag([0 0 -0.1]))
afc.addaniso('A')
afc.genmagstr('mode','direct','S',[0 0; 0 0; 1 -1]);
afc.field([0 0 7])
afcSpec = afc.spinwave({[0 0 0] [2 0 0] 101}, 'hermit',true);
figure; subplot(2,1,1)
sw_plotspec(afcSpec,'mode',4,'dE',0.2,'axLim',[0 3])
afcSpec = sw_egrid(sw_neutron(afcSpec),'Evect',linspace(0,6.5,500),'component','Sperp');
afcSpec = sw_omegasum(afcSpec,'zeroint',1e-6);
subplot(2,1,2)
sw_plotspec(afcSpec,'mode',2,'log',false,'axLim',[-4 10])
"""

if __name__ == "__main__":
    freeze_support()

    unit_cell = UnitCell(4,6,6)

    sites = [LatticeSite(0, 0, 0, 0,0,1, name="X"),
             LatticeSite(0.5,0,0, 0,0, -1, name="Y")]

    s = Structure(sites, unit_cell=unit_cell)

    exchanges = couplings(sites=sites,
                          unit_cell=unit_cell,
                          max_distance=1.1,
                          coupling_type=HeisenbergCoupling,
                          j=1,
                          direction_filter=filter([1,0,0], symmetric=True))

    anisotropies = axis_anisotropies(sites, -0.1)

    hamiltonian = Hamiltonian(s, exchanges, anisotropies)

    hamiltonian.print_summary()

    path = Path([[0,0,0], [2,0,0]])
    hamiltonian.plot(path, field=[0,0,7])
