""" Graphics representing a lattice site """
from pyspinw.gui.crystalviewer.GL.models import SolidModel
from pyspinw.site import LatticeSite

from pyspinw.gui.crystalviewer.GL.color import ColorSpecification, uniform_coloring
from pyspinw.gui.crystalviewer.GL.renderable import Renderable
from pyspinw.gui.crystalviewer.GL.transforms import Translation, Scaling
from pyspinw.gui.crystalviewer.arrowgraphics import ArrowGraphics
from pyspinw.gui.decorated import DecoratedSite



class SiteGraphics(Renderable):
    """ Renderable object representing lattice sites"""

    # Not selected or hover, red for now
    base_arrow = ArrowGraphics(uniform_coloring(1.0, 0.0, 0.0))
    group_base_arrow = ArrowGraphics(uniform_coloring(1.0, 0.3, 0.3))

    hover_arrow = ArrowGraphics(uniform_coloring(1.0, 0.6, 0.6))
    # Hover makes things brighter, hover and group is not a visually distinct case

    selected_arrow = ArrowGraphics(uniform_coloring(0.0, 1.0, 0.0))
    group_selected_arrow = ArrowGraphics(uniform_coloring(0.0, 1.0, 0.3))

    selected_hover_arrow = ArrowGraphics(uniform_coloring(0.6, 1.0, 0.6))
    # Again, no need for the group + hover

    # Group, Selected, Hover
    arrow_lookup = {
        (False, False, False): base_arrow,
        (False, False, True): hover_arrow,
        (False, True, False): selected_arrow,
        (False, True, True): selected_hover_arrow,
        (True, False, False): group_base_arrow,
        (True, False, True): hover_arrow,
        (True, True, False): group_selected_arrow,
        (True, True, True): selected_hover_arrow
    }


    def __init__(self, sites: list[DecoratedSite]):

        self._sites = sites
        self._scaling = 0.1
        self.graphics_objects = []

    @property
    def sites(self) -> list[DecoratedSite]:
        """ Get the current (decorated) sites"""
        return self._sites

    @sites.setter
    def sites(self, sites: list[DecoratedSite]):
        """ Set the sites to display """
        self._sites = sites
        self._update_site_graphics()

    def _update_site_graphics(self):
        # Build the tree
        new_graphics = []
        for site in self._sites:
            graphics = self.arrow_lookup[(site.flags.implied_hover, site.flags.selected, site.flags.hover)]

            new_graphics.append(Translation(
                site.position[0],
                site.position[1],
                site.position[2],
                Scaling(
                    self._scaling,
                    self._scaling,
                    self._scaling,
                graphics)))

        self.graphics_objects = new_graphics

    @property
    def scaling(self):
        """ Get the current arrow scaling"""
        return self._scaling

    @scaling.setter
    def scaling(self, scaling: float):
        """ Set the arrow scaling"""
        self._scaling = scaling
        self._update_site_graphics()


    def render_solid(self):
        """Renderable: how to render as a solid"""
        for graphics in self.graphics_objects:
            graphics.render_solid()

    def render_wireframe(self):
        """Renderable: how to render as a wireframe - does nothing right now"""
        return

