""" Kagome Ferromagnet example """

from multiprocessing.spawn import freeze_support

from pyspinw.coupling import HeisenbergCoupling
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.interface import spacegroup, couplings, filter
from pyspinw.path import Path, Path1D
from pyspinw.site import LatticeSite
from pyspinw.sample import Powder
from pyspinw.symmetry.supercell import TrivialSupercell
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structures import Structure
import sys

from pyspinw.debug_plot import debug_plot

"""
FMkagome = spinw;
FMkagome.genlattice('lat_const',[6 6 5],'angled',[90 90 120],'spgr','P -3')
FMkagome.addatom('r', [1/2 0 0], 'S', 1, 'label','MCu1','color','r')
FMkagome.gencoupling('maxDistance',4)
FMkagome.addmatrix('label','J1','value',-1,'color','orange');
FMkagome.addcoupling('mat','J1','bond',1);
FMkagome.genmagstr('mode','helical','k',[0 0 0],'n',[0 1 0],'S',[0 1 0]')
fmkSpec = FMkagome.spinwave({[-1/2 0 0] [0 0 0] [1/2 1/2 0] 100},'hermit',false);
fmkSpec = sw_egrid(sw_neutron(fmkSpec), 'Evect',linspace(0,6.5,100),'component','Sperp');
figure; sw_plotspec(fmkSpec,'mode',1,'colorbar',false,'axLim',[0 8])
fmkPow = FMkagome.powspec(linspace(0,2.5,100),'Evect',linspace(0,7,250),'nRand',1000,'hermit',false);
figure; sw_plotspec(fmkPow,'colorbar',true,'axLim',[0 0.05])
"""

if __name__ == "__main__":
    """Reproduces Tutorial 5: https://spinw.org/tutorials/05tutorial"""
    freeze_support()

    use_rust = "py" not in sys.argv[1] if len(sys.argv) > 1 else True

    unit_cell = UnitCell(6,6,5, gamma=120)

    x = LatticeSite(0.5, 0, 0, 0, 1, 0, name="X")
    y = LatticeSite(0, 0.5, 0, 0, 1, 0, name="Y")
    z = LatticeSite(0.5, 0.5, 0, 0, 1, 0, name="Z")

    sites = [x, y, z]

    s = Structure(sites, unit_cell=unit_cell, supercell=TrivialSupercell(scaling=(1,1,1)))


    exchanges = couplings(sites=[x, y, z],
                           unit_cell=unit_cell,
                           max_distance=4.,
                           coupling_type=HeisenbergCoupling,
                           j=-1)

    debug_plot(s, exchanges, show=False)

    hamiltonian = Hamiltonian(s, exchanges)

    hamiltonian.print_summary()

    path = Path([[-0.5,0,0], [0,0,0], [0.5,0.5,0]])
    hamiltonian.energy_plot(path, show=False, use_rust=use_rust)

    sample = Powder(hamiltonian)
    path1D = Path1D(0.1, 2.5, n_points=100)
    sample.show_spectrum(path1D, n_energy_bins=250)
