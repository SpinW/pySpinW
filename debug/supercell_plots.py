import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.lines import drawStyles
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

from pyspinw.site import LatticeSite
from pyspinw.symmetry.supercell import Supercell, TransformationSupercell, PropagationVector, RotationTransform, \
    CommensuratePropagationVector


def draw_supercell(figure: Figure, supercell: Supercell, sites: list[LatticeSite], vector_length: float=1.0, do_equal=False):
    ax = figure.add_subplot(projection='3d')

    xs = []
    ys = []
    zs = []
    mxs = []
    mys = []
    mzs = []
    for cell in supercell.cells():
        for site in sites:

            this_moment = supercell.moment(site, cell)

            xs.append(site.i + cell.i)
            ys.append(site.j + cell.j)
            zs.append(site.k + cell.k)

            mxs.append(this_moment[0])
            mys.append(this_moment[1])
            mzs.append(this_moment[2])

    ax.quiver(xs,ys,zs, mxs, mys, mzs, normalize=False, length=vector_length)
    if do_equal:
        ax.set_aspect("equal")


if __name__ == "__main__":
    fig = plt.figure("Test 1")

    supercell = TransformationSupercell([
        (CommensuratePropagationVector(1/3, 0, 0), RotationTransform([0,0,1])),
         (CommensuratePropagationVector(0, 1/3, 0) , RotationTransform([0,0,1])) ])

    sites = [
        LatticeSite(0.5, 0.5, 0.5, 1, 0, 0)
    ]

    draw_supercell(fig, supercell, sites, vector_length=0.3, do_equal=True)

    fig = plt.figure("Test 2")

    supercell = TransformationSupercell([
        (CommensuratePropagationVector(1 / 3, 0, 0), RotationTransform([0, 0, 1])),
        (CommensuratePropagationVector(0, 1 / 4, 0), RotationTransform([0, 0, 1])),
        (CommensuratePropagationVector(0, 0, 1/20), RotationTransform([0, 0, 1]))])

    sites = [
        LatticeSite(0.5, 0.5, 0.5, 1, 0, 0)
    ]

    draw_supercell(fig, supercell, sites, vector_length=0.3)

    plt.show()