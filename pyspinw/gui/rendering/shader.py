from abc import ABC, abstractmethod

from pyspinw.gui.rendering.load_shaders import load_shaders

from OpenGL.GL import *

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

