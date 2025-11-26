""" Temporary plot of sites and couplings """
import matplotlib.pyplot as plt

from pyspinw.coupling import Coupling
from pyspinw.site import LatticeSite
from pyspinw.structures import Structure


def debug_plot(structure: Structure, couplings: list[Coupling], show=True):
    fig = plt.figure()

    ax = fig.add_subplot(projection='3d')

    x,y,z = structure.matplotlib_site_data()

    ax.scatter(x,y,z)

    for coupling in couplings:
        pass

    if show:
        plt.show()