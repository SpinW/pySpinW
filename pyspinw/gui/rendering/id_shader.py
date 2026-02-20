""" Integer shader used for rendering object IDs as part of graphical selection """

import numpy as np
from OpenGL.GL import glGetUniformLocation, glUniformMatrix4fv, glUniform1ui, GL_TRUE

from pyspinw.gui.camera import Camera
from pyspinw.gui.rendering.shader import Shader


class IDShader(Shader):
    """ Integer shader used for rendering object IDs as part of graphical selection """

    vertex_shader = "simple_vertex"
    fragment_shader = "id_fragment"

    def __init__(self):
        super().__init__()

        self.model_matrix = np.eye(4, dtype=np.float32)
        self.projection_view = np.eye(4, dtype=np.float32)

        self.id_value = 0

        self._model_loc = glGetUniformLocation(self.shader_program, "model")
        self._projection_view_loc = glGetUniformLocation(self.shader_program, "projectionView")

        self._id_loc = glGetUniformLocation(self.shader_program, "objectID")



    @property
    def camera(self):
        """ Get the current camera"""
        return self._camera

    @camera.setter
    def camera(self, camera: Camera):
        """ Set the camera and associated matrices"""
        self._camera = camera

        view = self.camera.view_matrix()
        proj = self.camera.perspective_matrix(0.01, 100)
        self.projection_view = proj @ view

    def _set_uniforms(self):
        """ Set the shader uniforms for render"""
        glUniformMatrix4fv(self._model_loc, 1, GL_TRUE, np.array(self.model_matrix, dtype=np.float32))
        glUniformMatrix4fv(self._projection_view_loc, 1, GL_TRUE, np.array(self.projection_view, dtype=np.float32))

        glUniform1ui(self._id_loc, self.id_value)
