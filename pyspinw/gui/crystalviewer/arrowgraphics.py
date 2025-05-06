from pyspinw.gui.crystalviewer.GL.color import ColorSpecification
from pyspinw.gui.crystalviewer.GL.cone import Cone
from pyspinw.gui.crystalviewer.GL.cylinder import Cylinder
from pyspinw.gui.crystalviewer.GL.renderable import Renderable
from pyspinw.gui.crystalviewer.GL.transforms import Translation, Scaling


class ArrowGraphics(Renderable):
    def __init__(self,
                 color: ColorSpecification,
                 length: float = 1.0,
                 thickness: float = 1.0):

        cylinder_width_scale = 0.1 * thickness
        cone_scale = cylinder_width_scale * 3

        self.tree = Scaling(length, length, length,
            Translation(0,0,1,
                        Scaling(cone_scale, cone_scale, cone_scale,
                            Cone(colors=color))),
            Scaling(cylinder_width_scale, cylinder_width_scale, 1, Cylinder(colors=color)))

    def render_solid(self):
        self.tree.render_solid()