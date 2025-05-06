from typing import Callable

import numpy as np
from PySide6 import QtWidgets

from pyspinw.gui.crystalviewer.GL.color import uniform_coloring
from pyspinw.gui.crystalviewer.GL.scene import Scene
from pyspinw.gui.crystalviewer.arrowgraphics import ArrowGraphics
from pyspinw.gui.crystalviewer.unitcellgraphics import UnitCellGraphics
from pyspinw.unitcell import UnitCell

class CrystalViewer(Scene):
    def __init__(self, unit_cell: UnitCell,
                 parent=None,
                 on_key: Callable[[int], None] = lambda x: None):

        super().__init__(parent, on_key=on_key)

        self.unit_cell_graphics = UnitCellGraphics(unit_cell, uniform_coloring(1,1,1))

        self.add(self.unit_cell_graphics)

        # Set the scale and position of things

        self.view_centre = unit_cell.centre

        self.min_distance = 0.01
        self.max_distance = 250

        self.view_distance = 3*np.linalg.norm(unit_cell.centre)


def main():
    """ Show a demo of the opengl.rst window """
    import os

    os.environ["QT_ENABLE_HIGHDPI_SCALING"] = "1"
    app = QtWidgets.QApplication([])

    mainWindow = QtWidgets.QMainWindow()


    unit_cell = UnitCell(3,4,5, 60, 45)
    viewer = CrystalViewer(unit_cell)

    viewer.add(ArrowGraphics(uniform_coloring(1,0,0)))

    mainWindow.setCentralWidget(viewer)
    mainWindow.show()

    mainWindow.resize(600, 600)
    app.exec_()


if __name__ == "__main__":
    main()