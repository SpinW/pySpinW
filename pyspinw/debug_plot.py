"""Temporary plot of sites and couplings."""
import matplotlib.pyplot as plt

from pyspinw.exchange import Exchange
from pyspinw.site import Site
from pyspinw.structure import Structure


def debug_plot(structure: Structure, couplings: list[Exchange], show=True):
    """Plotting.

    Parameters
    ----------
    structure : Structure
        Crystal structure containing sites to plot.
    couplings : list of Exchange
        List of exchange couplings to plot.
    show : bool, optional
        If ``True``, call `` matplotlib.pyplot.show()``. Default is ``True``.
    """
    fig = plt.figure()

    ax = fig.add_subplot(projection='3d')

    x,y,z = structure.matplotlib_site_data()

    ax.scatter(x,y,z)

    for coupling in couplings:
        pass

    ax.axis("equal")

    if show:
        plt.show()
