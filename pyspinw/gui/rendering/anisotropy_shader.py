""" Shader for most of the objects in the viewer"""

import numpy as np
from OpenGL.GL import *

from pyspinw.gui.camera import Camera
from pyspinw.gui.rendering.shader import Shader


class AnisotropyShader(Shader):
    """ Shader for most of the objects in the viewer"""

    vertex_shader = "object_vertex"
    fragment_shader = "fresnel_fragment"

    def __init__(self):
        super().__init__()

        self._camera = None

        self._model_loc = glGetUniformLocation(self.shader_program, "model")
        self._projection_view_loc = glGetUniformLocation(self.shader_program, "projectionView")

        self._view_pos_loc = glGetUniformLocation(self.shader_program, "viewPos")

        self._object_color_loc = glGetUniformLocation(self.shader_program, "objectColor")

        self.model_matrix = np.eye(4, dtype=np.float32)
        self.projection_view = np.eye(4, dtype=np.float32)
        self.view = np.eye(4, dtype=np.float32)

        self.view_position = 0, 0, 0

        self.object_color = 1, 1, 1


    @property
    def camera(self):
        """ Get the current camera"""
        return self._camera

    @camera.setter
    def camera(self, camera: Camera):
        """ Set the camera and associated matrices """
        self._camera = camera

        self.view = self.camera.view_matrix()
        proj = self.camera.projection_matrix(0.01, 100)
        self.projection_view = proj @ self.view

        self.view_position = camera.position

    def _set_uniforms(self):

        glUniformMatrix4fv(self._model_loc, 1, GL_TRUE, np.array(self.model_matrix, dtype=np.float32))
        glUniformMatrix4fv(self._projection_view_loc, 1, GL_TRUE, np.array(self.projection_view, dtype=np.float32))

        glUniform3f(self._view_pos_loc, *self.view_position)

        glUniform3f(self._object_color_loc, *self.object_color)


