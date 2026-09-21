""" Twinned Ferromagnetic chain example """
import numpy as np
from pyspinw import *

unit_cell = UnitCell(1,1,1)

only_site = LatticeSite(0, 0, 0, 0,0,1, name="X")

s = Structure([only_site], unit_cell=unit_cell)

exchanges = generate_exchanges(sites=[only_site],
                               unit_cell=unit_cell,
                               max_distance=1.1,
                               exchange_type=HeisenbergExchange,
                               j=-1,
                               direction_filter=filter([1,0,0]))

hamiltonian = Hamiltonian(s, exchanges)

path = Path([[0,0,0], [1,0,0]])
y_path = Path([[0,0,0],[0,1,0]])


single_crystal = SingleCrystal(hamiltonian)
single_crystal.spaghetti_plot(path, dE=0.4)

single_crystal.spaghetti_plot(y_path, dE=0.4, evect=np.linspace(0,4,100))

twin = Twin(hamiltonian, np.array([[0,1,0],[-1,0,0],[0,0,1]], dtype=float), 0.4)

twin.spaghetti_plot(path, dE=0.4)

