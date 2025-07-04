""" Transformations of objects """

import logging
from typing import List

import numpy as np

from OpenGL.GL import *

from pyspinw.gui.crystalviewer.GL.renderable import Renderable


logger = logging.getLogger("GL.transforms")

class SceneGraphNode(Renderable):
    """General transform - also doubles as a scene graph node

    For the sake of speed, the transformation matrix shape is not checked.
    It should be a 4x4 transformation matrix
    """

    def __init__(self, *children: Renderable):
        super().__init__()
        self.children: List[Renderable] = list(children)
        self.solid_render_enabled = True
        self.wireframe_render_enabled = True


    def add_child(self, child: Renderable):
        """ Add a renderable object to this scene graph node"""
        self.children.append(child)

    def apply(self):
        """ GL operations needed to apply any transformations associated with this node """

    def render_solid(self):
        """Renderable: solid rendering implementation"""
        if self.solid_render_enabled:

            # Apply transform
            glPushMatrix()
            self.apply()

            # Check stack
            if glGetIntegerv(GL_MODELVIEW_STACK_DEPTH) >= 16:
                logger.info(f"GL Stack size utilisation {glGetIntegerv(GL_MODELVIEW_STACK_DEPTH)}" +
                            "the limit could be as low as 16 on some machines")

            # Render children
            for child in self.children:
                child.render_solid()

            # Unapply
            glPopMatrix()


    def render_wireframe(self):
        """Renderable: wireframe rendering implementation"""
        if self.wireframe_render_enabled:
            # Apply transform
            glPushMatrix()
            self.apply()

            # Check stack
            if glGetIntegerv(GL_MODELVIEW_STACK_DEPTH) >= 16:
                logger.info(f"GL Stack size utilisation {glGetIntegerv(GL_MODELVIEW_STACK_DEPTH)}, "
                            "the limit could be as low as 16 on some machines")


            # Render children
            for child in self.children:
                child.render_wireframe()

            # unapply transform
            glPopMatrix()


class Rotation(SceneGraphNode):
    """ Rotation of a renderable object"""

    def __init__(self, angle, x, y, z, *children: Renderable):
        """Rotate the children of this node

        :param angle: angle of rotation in degrees
        :param axis: axis for rotation
        """
        super().__init__(*children)
        self.angle = angle
        self.x = x
        self.y = y
        self.z = z

    def apply(self):
        """Transform implementation: gl operation that can be popped"""
        glRotate(self.angle, self.x, self.y, self.z)


class Translation(SceneGraphNode):
    """ Translation of a renderable object"""

    def __init__(self, x: float, y: float, z: float, *children: Renderable):
        """Translate the children of this node

        :param x: x translation
        :param y: y translation
        :param z: z translation
        """
        super().__init__(*children)
        self.x = x
        self.y = y
        self.z = z

    def apply(self):
        """Transform implementation: gl operation that can be popped"""
        glTranslate(self.x, self.y, self.z)


class Scaling(SceneGraphNode):
    """ Scaling of a renderable object"""

    def __init__(self, x: float, y: float, z: float, *children: Renderable):
        """Scale the children of this node

        :param x: x scale
        :param y: y scale
        :param z: z scale
        """
        super().__init__(*children)
        self.x = x
        self.y = y
        self.z = z

    def apply(self):
        """Transform implementation: gl operation that can be popped"""
        glScale(self.x, self.y, self.z)

class MatrixTransform(SceneGraphNode):
    """ General matrix transform of a renderable object"""

    def __init__(self, matrix: np.ndarray, *children: Renderable):
        """Apply a 4x4 transformation matrix to the children of this node

        :param matrix: a 4x4 transformation matrix

        """
        super().__init__(*children)

        self.matrix = matrix

    def apply(self):
        """Transform implementation: gl operation that can be popped"""
        glMultMatrixd(self.matrix)
