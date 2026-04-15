""" Problem 5: Samples"""

# Part 1: Make a twin

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

twin = Twin(hamiltonian, np.array([[0,1,0],[-1,0,0],[0,0,1]], dtype=float), 0.4)

twin.spaghetti_plot(path, dE=0.4)


#
# Part 2: 'sall good, man!
#