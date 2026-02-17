import numpy as np

from pyspinw.gui.camera import Camera
from pyspinw.gui.rendering.shader import Shader

from pyspinw.gui.rendering.renderable_objects import SelectionMode

from OpenGL.GL import *

class CellShader(Shader):
    vertex_shader = "simple_vertex"
    fragment_shader = "cell_fragment"

    def __init__(self):
        super().__init__()

        self.mode = SelectionMode.NOT_SELECTED

        self._model_loc = glGetUniformLocation(self.shader_program, "model")
        self._projection_view_loc = glGetUniformLocation(self.shader_program, "projectionView")

        self.model_matrix = np.eye(4, dtype=np.float32)
        self.projection_view = np.eye(4, dtype=np.float32)

    @property
    def camera(self):
        return self._camera

    @camera.setter
    def camera(self, camera: Camera):
        self._camera = camera

        view = self.camera.view_matrix()
        proj = self.camera.perspective_matrix(0.01, 100)
        self.projection_view = proj @ view

    def _set_uniforms(self):

        glUniformMatrix4fv(self._model_loc, 1, GL_TRUE, np.array(self.model_matrix, dtype=np.float32))
        glUniformMatrix4fv(self._projection_view_loc, 1, GL_TRUE, np.array(self.projection_view, dtype=np.float32))
