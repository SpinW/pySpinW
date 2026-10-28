"""
Functions to switch a default color scheme of plots to dark mode
"""

from matplotlib.spines import Spine
from matplotlib.patches import Polygon, Rectangle
from matplotlib.text import Text
from pyspinw.gui.displayoptions import DisplayOptions
from pyspinw.gui.viewer import snapshot

def _set_white_recursive(artist):
    if isinstance(artist, (Spine, Polygon, Rectangle)):
        artist.set_edgecolor('white')
    elif isinstance(artist, (Text,)):
        artist.set_color('white')
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
