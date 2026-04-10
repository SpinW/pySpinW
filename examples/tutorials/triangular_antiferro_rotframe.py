""" Triangular Antiferromagnet example

Reproduces Tutorial 12: https://spinw.org/tutorials/12tutorial"""

from pyspinw import *
import matplotlib.pyplot as plt
from pyspinw.legacy.genmagstr import genmagstr

unit_cell = UnitCell(3, 3, 8, gamma=120)

sites = generate_helical_structure(unit_cell, positions=[[0, 0, 0], [0, 0, 0.5]], spins=[[0, 1, 0], [0, 1, 0]],
                                   magnitudes=[3./2, 3./2], propagation_vector=[1./3, 1./3, 0], perpendicular=[0, 0, 1])

exchanges = generate_exchanges(sites=sites, bond=1, exchange_type=HeisenbergExchange, j=1) \
          + generate_exchanges(sites=sites, bond=2, exchange_type=HeisenbergExchange, j=-0.1)

anisotropies = axis_anisotropies(sites, 0.2)
hamiltonian = Hamiltonian(sites, exchanges, anisotropies)

hamiltonian.print_summary()

path = Path([[0,0,0], [1,1,0]], n_points_per_segment=401)

fig = hamiltonian.spaghetti_plot(path, show=False, use_rotating=True)
fig.axes[0].set_ylim(0, 6)
plt.show()
