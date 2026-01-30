from abc import ABC, abstractmethod

import numpy as np
from OpenGL.GL.shaders import GL_FALSE
from OpenGL.raw.GL.VERSION.GL_2_0 import (glUseProgram, glGetUniformLocation, glUniform3f,
                                          glUniformMatrix4fv, glUniform3fv, glUniform1f, GL_TRUE)

from pyspinw.gui.camera import Camera
from pyspinw.gui.rendering.load_shaders import load_shaders
from pyspinw.gui.rendering.renderable_objects import SelectionMode


class Shader(ABC):
    """ Base class for shader objects

    These objects must be initialised in the initialiseGL section of an OpenGLWidget, or, at least, after open gl
    has been started up.
    """

    vertex_shader = "default_vertex"
    fragment_shader = "default_fragment"

    def __init__(self):
        self.shader_program = load_shaders(fragment_filename=self.fragment_shader, vertex_filename=self.vertex_shader)

    @abstractmethod
    def _set_uniforms(self):
        """ Set up shader variables"""

    def use(self):
        """ Call to set this shader for use"""

        glUseProgram(self.shader_program)
        self._set_uniforms()


class SelectionShader(Shader):
    """ Used for rendering the selection highlighting """

    vertex_shader = "selection_vertex"
    fragment_shader = "selection_fragment"

    def __init__(self):
        super().__init__()

        self.mode = SelectionMode.NOT_SELECTED
        self._selection_color_loc = glGetUniformLocation(self.shader_program, "selectionColor")

    def _set_uniforms(self):
        """ Set the shader uniforms for render"""
        match self.mode:
            case SelectionMode.NOT_SELECTED:
                glUniform3f(self._selection_color_loc, 1.0, 0.2, 0.2)
            case SelectionMode.SELECTED:
                glUniform3f(self._selection_color_loc, 1.0, 0.6, 0.1)
            case SelectionMode.HOVER:
                glUniform3f(self._selection_color_loc, 1.0, 1.0, 1.0)
            case SelectionMode.SECONDARY:
                glUniform3f(self._selection_color_loc, 0.8, 0.8, 0.1)


class ObjectShader(Shader):

    vertex_shader = "object_vertex"
    fragment_shader = "object_fragment"

    def __init__(self):
        super().__init__()

        self._camera = None

        self._model_loc = glGetUniformLocation(self.shader_program, "model")
        self._projection_view_loc = glGetUniformLocation(self.shader_program, "projectionView")

        self._light_pos_loc = glGetUniformLocation(self.shader_program, "lightPos")
        self._view_pos_loc = glGetUniformLocation(self.shader_program, "viewPos")

        self._light_color_loc = glGetUniformLocation(self.shader_program, "lightColor")
        self._object_color_loc = glGetUniformLocation(self.shader_program, "objectColor")

        self._specular_strength_loc = glGetUniformLocation(self.shader_program, "specularStrength")
        self._ambient_strength_loc = glGetUniformLocation(self.shader_program, "ambientStrength")

        self.model_matrix = np.eye(4, dtype=np.float32)
        self.projection_view = np.eye(4, dtype=np.float32)

        self.light_position = -20, 0, 0
        self.view_position = 0, 0, 0

        self.light_color = 1, 1, 1
        self.object_color = 1, 1, 1

        self.specular_strength = 0.3
        self.ambient_strength = 0.2

    @property
    def camera(self):
        return self._camera

    @camera.setter
    def camera(self, camera: Camera):
        self._camera = camera

        view = self.camera.view_matrix()
        proj = self.camera.perspective_matrix(0.01, 100)
        self.projection_view = proj @ view

        self.view_position = camera.position

    def _set_uniforms(self):
        glUniformMatrix4fv(self._model_loc, 1, GL_FALSE, np.array(self.model_matrix, dtype=np.float32))
        glUniformMatrix4fv(self._projection_view_loc, 1, GL_TRUE, np.array(self.projection_view, dtype=np.float32))

        glUniform3f(self._light_pos_loc, *self.light_position)
        glUniform3f(self._view_pos_loc, *self.view_position)

        glUniform3f(self._light_color_loc, *self.light_color)
        glUniform3f(self._object_color_loc, *self.object_color)

        glUniform1f(self._specular_strength_loc, self.specular_strength)
        glUniform1f(self._ambient_strength_loc, self.ambient_strength)
