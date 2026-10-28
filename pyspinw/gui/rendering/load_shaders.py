"""Helper methods for loading and compiling shaders."""

import logging
import os.path
from importlib import resources

from OpenGL.GL import (glCreateShader, glShaderSource, glCompileShader, glGetShaderiv, glGetShaderInfoLog,
                       glCreateProgram, glAttachShader, glLinkProgram, glGetProgramiv, glGetProgramInfoLog,
                       glDeleteShader, GL_COMPILE_STATUS, GL_VERTEX_SHADER, GL_FRAGMENT_SHADER, GL_LINK_STATUS)



logger = logging.Logger("pyspinw.gui.rendering.shaders")

def _compile_shader(source, shader_type):
    """Compile a shader.

    Parameters
    ----------
    source : str
        GLSL shader source code.
    shader_type : int
        OpenGL shader type constant (e.g., ``GL_VERTEX_SHADER`` or ``GL_FRAGMENT_SHADER``).

    Returns
    -------
    int
        Compiled OpenGL shader handle.

    Raises
    ------
    RuntimeError
        If shader compilation fails.
    """
    shader = glCreateShader(shader_type)

    glShaderSource(shader, source)
    glCompileShader(shader)

    success = glGetShaderiv(shader, GL_COMPILE_STATUS)
    if not success:
        info_log = glGetShaderInfoLog(shader).decode()
        raise RuntimeError(f"Shader compile error: {info_log}")

    return shader

def _create_shader_program(vs_src, fs_src):
    """Create a program object with vertex and fragment shaders.

    Parameters
    ----------
    vs_src : str
        Vertex shader GLSL source code.
    fs_src : str
        Fragment shader GLSL source code.

    Returns
    -------
    int
        Linked OpenGL shader program handle.

    Raises
    ------
    RuntimeError
        If shader compilation or program linking fails.
    """
    vertex_shader = _compile_shader(vs_src, GL_VERTEX_SHADER)
    fragment_shader = _compile_shader(fs_src, GL_FRAGMENT_SHADER)

    program = glCreateProgram()

    glAttachShader(program, vertex_shader)
    glAttachShader(program, fragment_shader)
    glLinkProgram(program)

    success = glGetProgramiv(program, GL_LINK_STATUS)
    if not success:
        info_log = glGetProgramInfoLog(program).decode()
        raise RuntimeError(f"GL program link error: {info_log}")

    glDeleteShader(vertex_shader)
    glDeleteShader(fragment_shader)

    return program


def load_shaders(fragment_filename: str | None = None, vertex_filename: str | None = None):
    """Load shaders from file, or from the shaders directory.

    Parameters
    ----------
    fragment_filename : str or None, optional
        File path or shader resource name for the fragment shader. Default is ``None``.
    vertex_filename : str or None, optional
        File path or shader resource name for the vertex shader. Default is ``None``.

    Returns
    -------
    int
        Linked OpenGL shader program handle.
    """
    if fragment_filename is None:
        frag_code = resources.read_text("pyspinw.gui.rendering.shaders", "default_fragment.glsl")

    else:
        if os.path.exists(fragment_filename):
            with open(fragment_filename, 'r') as file:
                frag_code = "".join(file.readlines())
        else:
            frag_code = resources.read_text("pyspinw.gui.rendering.shaders", fragment_filename + ".glsl")

    if vertex_filename is None:
        vertex_code = resources.read_text("pyspinw.gui.rendering.shaders", "default_vertex.glsl")

    else:
        if os.path.exists(vertex_filename):
            with open(vertex_filename, 'r') as file:
                vertex_code = "".join(file.readlines())
        else:
            vertex_code = resources.read_text("pyspinw.gui.rendering.shaders", vertex_filename + ".glsl")

    return _create_shader_program(vertex_code, frag_code)
