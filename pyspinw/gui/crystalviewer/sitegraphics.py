""" Graphics representing a lattice site """
from pyspinw.gui.crystalviewer.GL.renderable import Renderable
from pyspinw.site import LatticeSite


class SpinSite(Renderable):
    def __init__(self, site: LatticeSite):
        self.site = site

    def render_solid(self):
        pass