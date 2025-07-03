from typing import Callable

import numpy as np
from PySide6 import QtWidgets

from pyspinw.gui.crystalviewer.GL.color import uniform_coloring
from pyspinw.gui.crystalviewer.GL.scene import Scene
from pyspinw.gui.crystalviewer.arrowgraphics import ArrowGraphics
from pyspinw.gui.crystalviewer.sitegraphics import SiteGraphics
from pyspinw.gui.crystalviewer.unitcellgraphics import UnitCellGraphics
from pyspinw.gui.decorated import DecoratedSite
from pyspinw.symmetry.unitcell import UnitCell

class CrystalViewer(Scene):
    def __init__(self,
                 parent=None,
                 on_key: Callable[[int], None] = lambda x: None):

        super().__init__(parent, on_key=on_key)

        # Unit cell

        self._unit_cell = None
        self._unit_cell_graphics = UnitCellGraphics(UnitCell(1,1,1,90,90,90), uniform_coloring(1, 1, 1))
        self._unit_cell_graphics.wireframe_render_enabled = False

        self.add(self._unit_cell_graphics)

        # Sites
        self._site_graphics = SiteGraphics(sites=[])
        self.add(self._site_graphics)

        # Set the scale and position of things

        self.min_distance = 0.01
        self.max_distance = 250


    @property
    def unit_cell(self) -> UnitCell | None:
        return self.unit_cell

    @unit_cell.setter
    def unit_cell(self, unit_cell: UnitCell | None):
        self._unit_cell = unit_cell

        if unit_cell is None:
            self._unit_cell_graphics.wireframe_render_enabled = False
        else:
            self._unit_cell_graphics.wireframe_render_enabled = True
            self._unit_cell_graphics.cell = unit_cell
            self.view_centre = unit_cell.centre
            self.view_distance = 3*np.linalg.norm(unit_cell.centre)

        self.update()

    @property
    def sites(self) -> list[DecoratedSite]:
        pass

    @sites.setter
    def sites(self, sites: list[DecoratedSite]):
        self._site_graphics.sites = sites

def main():
    """ Show a demo of the opengl window """
    import os

    os.environ["QT_ENABLE_HIGHDPI_SCALING"] = "1"
    app = QtWidgets.QApplication([])

    mainWindow = QtWidgets.QMainWindow()


    unit_cell = UnitCell(3,4,5, 60, 45)
    viewer = CrystalViewer()

    viewer.add(ArrowGraphics(uniform_coloring(1,0,0)))

    mainWindow.setCentralWidget(viewer)
    mainWindow.show()

    mainWindow.resize(600, 600)
    app.exec_()


if __name__ == "__main__":
    main()
