""" Graphical representation of the unit cell"""

import numpy as np


from pyspinw.gui.crystalviewer.GL.color import ColorSpecification
from pyspinw.gui.crystalviewer.GL.models import WireModel
from pyspinw.gui.crystalviewer.GL.cube import Cube

from pyspinw.symmetry.unitcell import UnitCell



class UnitCellGraphics(WireModel):
    """ Graphical representation of the unit cell"""

    def __init__(self,
                 unit_cell: UnitCell,
                 edge_colors: ColorSpecification | None):

        self._cell = unit_cell

        # Based on cube primitive geometry
        vertices = [unit_cell.fractional_to_cartesian(np.array([vert])+0.5)[0,:] for vert in Cube.cube_vertices]

        edges = Cube.cube_edges

        super().__init__(
            vertices=vertices,
            edges=edges,
            edge_colors=edge_colors)

    @property
    def cell(self):
        return self._cell

    @cell.setter
    def cell(self, unit_cell: UnitCell):
        self._cell = unit_cell
        self.vertices = [unit_cell.fractional_to_cartesian(np.array([vert])+0.5)[0,:] for vert in Cube.cube_vertices]
