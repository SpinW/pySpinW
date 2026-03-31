""" Triangular Antiferromagnet example

Reproduces Tutorial 12: https://spinw.org/tutorials/12tutorial"""

from pyspinw import *
import matplotlib.pyplot as plt
from pyspinw.legacy.genmagstr import genmagstr

unit_cell = UnitCell(3, 3, 8, gamma=120)

sites = [LatticeSite(0, 0, 0, 0, 1, 0, name="X"),
         LatticeSite(0, 0, 0.5, 0, 1, 0, name="Y")]

s = genmagstr(sites, unit_cell, magnitude=[3./2, 3./2],mode='helical', k=[1./3, 1./3, 0], n=[0, 0, 1])

exchanges = generate_exchanges(sites=sites, unit_cell=unit_cell, max_distance=3.1, coupling_type=HeisenbergCoupling, j=1) \
          + generate_exchanges(sites=sites, unit_cell=unit_cell, min_distance=3.1, max_distance=4.1, j=-0.1, coupling_type=HeisenbergCoupling)

anisotropies = axis_anisotropies(sites, 0.2)
hamiltonian = Hamiltonian(s, exchanges, anisotropies)

hamiltonian.print_summary()

path = Path([[0,0,0], [1,1,0]], n_points_per_segment=401)

fig = hamiltonian.spaghetti_plot(path, show=False, use_rotating=True)
fig.axes[0].set_ylim(0, 7)
fig.axes[1].set_ylim(0, 5)
plt.show()
