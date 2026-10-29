"""
Functions to switch a default color scheme of plots to dark mode
"""

from matplotlib.spines import Spine
from matplotlib.patches import Polygon, Rectangle
from matplotlib.text import Text
from matplotlib.lines import Line2D
from matplotlib.collections import QuadMesh
from pyspinw.gui.displayoptions import DisplayOptions
from pyspinw.gui.viewer import snapshot
import numpy as np


def _set_white_recursive(artist):
    print(type(artist), artist)
    if isinstance(artist, QuadMesh):
        colors = artist.get_cmap().colors
        # Check if colormap goes from light to dark (light==sqrt(3)==1.7)
        if np.linalg.norm(colors[0]) > 1.5 and np.linalg.norm(colors[-1]) < 0.1:
            artist.set_cmap(artist.get_cmap().reversed())
    if isinstance(artist, (Spine, Polygon, Rectangle)):
        artist.set_edgecolor('white')
    elif isinstance(artist, (Text,)):
        artist.set_color('white')
    elif isinstance(artist, (Line2D,)):
        color = artist.get_color()
        if color[0:3] == (0, 0, 0):
            artist.set_color((1, 1, 1, color[3]))
    elif hasattr(artist, 'get_children'):
        for child in artist.get_children():
            _set_white_recursive(child)


def savefig_dark(fig, outputfile):
    fig.set_facecolor('black')
    for ax in fig.axes:
        ax.set_facecolor('black')
        _set_white_recursive(ax)
    fig.savefig(outputfile)


def snapshot_light(*args, **kwargs):
    kwargs['display_options'] = DisplayOptions.LightMode()
    snapshot(*args, **kwargs)


def snapshot_dark(*args, **kwargs):
    kwargs['display_options'] = DisplayOptions.DarkMode()
    snapshot(*args, **kwargs)
