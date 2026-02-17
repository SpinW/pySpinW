import numpy as np
from OpenGL.GL import *

from pyspinw.gui.camera import Camera
from pyspinw.gui.rendering.renderable_objects import SelectionMode
from pyspinw.gui.rendering.shader import Shader


class SelectionShader(Shader):
    """ Used for rendering the selection highlighting """

    vertex_shader = "simple_vertex"
    fragment_shader = "selection_fragment"

    def __init__(self):
        super().__init__()

        self.mode = SelectionMode.NOT_SELECTED

        self.model_matrix = np.eye(4, dtype=np.float32)
        self.projection_view = np.eye(4, dtype=np.float32)

        # Uniform locations
        self._model_loc = glGetUniformLocation(self.shader_program, "model")
        self._projection_view_loc = glGetUniformLocation(self.shader_program, "projectionView")

        self._selection_color_loc = glGetUniformLocation(self.shader_program, "selectionColor")


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
        """ Set the shader uniforms for render"""

        glUniformMatrix4fv(self._model_loc, 1, GL_TRUE, np.array(self.model_matrix, dtype=np.float32))
        glUniformMatrix4fv(self._projection_view_loc, 1, GL_TRUE, np.array(self.projection_view, dtype=np.float32))

        match self.mode:
            case SelectionMode.NOT_SELECTED:
                glUniform3f(self._selection_color_loc, 1.0, 0.2, 0.2)
            case SelectionMode.SELECTED:
                glUniform3f(self._selection_color_loc, 1.0, 0.6, 0.1)
            case SelectionMode.HOVER:
                glUniform3f(self._selection_color_loc, 1.0, 1.0, 1.0)
            case SelectionMode.SECONDARY:
                glUniform3f(self._selection_color_loc, 0.8, 0.8, 0.1)

