"""Shader used for highlighting."""

import numpy as np
from OpenGL.GL import *

from pyspinw.gui.camera import Camera
from pyspinw.gui.rendering.selectionmode import SelectionMode
from pyspinw.gui.rendering.shader import Shader


class SelectionShader(Shader):
    """Used for rendering the selection highlighting.

    Attributes
    ----------
    vertex_shader : str
        Name of the vertex shader resource.
    fragment_shader : str
        Name of the fragment shader resource.
    mode : SelectionMode
        Current selection highlighting state.
    model_matrix : (4, 4) numpy.ndarray
        Model transformation matrix.
    projection_view : (4, 4) numpy.ndarray
        Combined projection-view matrix.
    selected_color : tuple of (float, float, float)
        RGB highlight color for selected objects.
    hover_color : tuple of (float, float, float)
        RGB highlight color for hovered objects.
    selected_hover_color : tuple of (float, float, float)
        RGB highlight color for objects that are both selected and hovered.
    """

    vertex_shader = "simple_vertex"
    fragment_shader = "emission_fragment"

    def __init__(self):
        super().__init__()

        self.mode = SelectionMode.NOT_SELECTED

        self.model_matrix = np.eye(4, dtype=np.float32)
        self.projection_view = np.eye(4, dtype=np.float32)

        # Colours for selections
        self.selected_color = 1.0, 0.6, 0.1
        self.hover_color = 1.0, 1.0, 1.0
        self.selected_hover_color = 1.0, 0.8, 0.1

        # Uniform locations
        self._model_loc = glGetUniformLocation(self.shader_program, "model")
        self._projection_view_loc = glGetUniformLocation(self.shader_program, "projectionView")

        self._selection_color_loc = glGetUniformLocation(self.shader_program, "objectColor")



    @property
    def camera(self):
        """Current camera used by the shader."""
        return self._camera

    @camera.setter
    def camera(self, camera: Camera):
        """Set the camera associated matrices."""
        self._camera = camera

        view = self.camera.view_matrix()
        proj = self.camera.projection_matrix(0.01, 100)
        self.projection_view = proj @ view

    def _set_uniforms(self):
        """Set the shader uniforms for render."""
        glUniformMatrix4fv(self._model_loc, 1, GL_TRUE, np.array(self.model_matrix, dtype=np.float32))
        glUniformMatrix4fv(self._projection_view_loc, 1, GL_TRUE, np.array(self.projection_view, dtype=np.float32))

        match self.mode:
            case SelectionMode.NOT_SELECTED:
                glUniform3f(self._selection_color_loc, 1.0, 1.0, 1.0)
            case SelectionMode.SELECTED:
                glUniform3f(self._selection_color_loc, *self.selected_color)
            case SelectionMode.HOVER:
                glUniform3f(self._selection_color_loc, *self.hover_color)
            case SelectionMode.SELECTED_HOVER:
                glUniform3f(self._selection_color_loc, *self.selected_hover_color)

