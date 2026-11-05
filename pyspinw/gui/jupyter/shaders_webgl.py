"""WebGL 2.0 (GLSL ES 3.00) shaders and shader program managers for Jupyter notebooks."""

from __future__ import annotations

from typing import TYPE_CHECKING
import numpy as np

if TYPE_CHECKING:
    from pyspinw.gui.camera import Camera

# ---------------------------------------------------------------------------
# GLSL ES 3.00 Shader Sources
# ---------------------------------------------------------------------------

OBJECT_VERTEX_SHADER = """#version 300 es
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;

uniform mat4 model;
uniform mat4 projectionView;

out vec3 FragPos;
out vec3 Normal;

void main()
{
    FragPos = vec3(model * vec4(aPos, 1.0));
    Normal = mat3(transpose(inverse(model))) * aNormal;
    gl_Position = projectionView * vec4(FragPos, 1.0);
}
"""

OBJECT_FRAGMENT_SHADER = """#version 300 es
precision highp float;

in vec3 FragPos;
in vec3 Normal;

out vec4 FragColor;

uniform vec3 lightPos;
uniform vec3 viewPos;
uniform vec3 lightColor;
uniform vec3 objectColor;
uniform float ambientStrength;
uniform float specularStrength;

void main()
{
    vec3 ambient = ambientStrength * lightColor;

    vec3 norm = normalize(Normal);
    vec3 lightDir = normalize(lightPos - FragPos);
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = diff * lightColor;

    vec3 viewDir = normalize(viewPos - FragPos);
    vec3 halfwayDir = normalize(lightDir + viewDir);
    float spec = pow(max(dot(norm, halfwayDir), 0.0), 32.0);
    vec3 specular = specularStrength * spec * lightColor;

    vec3 result = (ambient + diffuse) * objectColor + specular;
    FragColor = vec4(result, 1.0);
}
"""

SIMPLE_VERTEX_SHADER = """#version 300 es
layout (location = 0) in vec3 aPos;

uniform mat4 model;
uniform mat4 projectionView;

void main()
{
    gl_Position = projectionView * model * vec4(aPos, 1.0);
}
"""

CELL_FRAGMENT_SHADER = """#version 300 es
precision highp float;

uniform vec3 cellColor;
out vec4 FragColor;

void main()
{
    FragColor = vec4(cellColor, 1.0);
}
"""

FRESNEL_FRAGMENT_SHADER = """#version 300 es
precision highp float;

in vec3 FragPos;
in vec3 Normal;

out vec4 FragColor;

uniform vec3 viewPos;
uniform vec3 objectColor;

void main()
{
    vec3 N = normalize(Normal);
    vec3 V = normalize(viewPos - FragPos);
    float alpha = pow(1.0 - max(dot(N, V), 0.0), 2.0);
    FragColor = vec4(objectColor, alpha);
}
"""

EMISSION_FRAGMENT_SHADER = """#version 300 es
precision highp float;

uniform vec3 objectColor;
out vec4 FragColor;

void main()
{
    FragColor = vec4(objectColor, 1.0);
}
"""


def compile_webgl_program(viewer, vs_src: str, fs_src: str):
    """Compile and link a WebGL2 program on an ipywebgl.GLViewer.

    Parameters
    ----------
    viewer : ipywebgl.GLViewer
        The WebGL viewer widget.
    vs_src : str
        Vertex shader GLSL ES 3.00 source code.
    fs_src : str
        Fragment shader GLSL ES 3.00 source code.

    Returns
    -------
    program
        The linked WebGL program handle/object.
    """
    vs = viewer.create_shader("VERTEX_SHADER")
    viewer.shader_source(vs, vs_src)
    viewer.compile_shader(vs)

    fs = viewer.create_shader("FRAGMENT_SHADER")
    viewer.shader_source(fs, fs_src)
    viewer.compile_shader(fs)

    program = viewer.create_program()
    viewer.attach_shader(program, vs)
    viewer.attach_shader(program, fs)
    viewer.link_program(program)

    return program


def _mat_to_col_major(matrix: np.ndarray) -> np.ndarray:
    """Ensure matrix is flattened in column-major order for WebGL."""
    return np.ascontiguousarray(matrix.T, dtype=np.float32).flatten()


# ---------------------------------------------------------------------------
# High-Level Shader Program Wrappers
# ---------------------------------------------------------------------------

class ObjectShaderWebGL:
    """WebGL2 shader for atoms, spins, and exchange bonds."""

    def __init__(self, viewer):
        self.viewer = viewer
        self.program = compile_webgl_program(viewer, OBJECT_VERTEX_SHADER, OBJECT_FRAGMENT_SHADER)

        self._model_loc = viewer.get_uniform_location(self.program, "model")
        self._projection_view_loc = viewer.get_uniform_location(self.program, "projectionView")
        self._light_pos_loc = viewer.get_uniform_location(self.program, "lightPos")
        self._view_pos_loc = viewer.get_uniform_location(self.program, "viewPos")
        self._light_color_loc = viewer.get_uniform_location(self.program, "lightColor")
        self._object_color_loc = viewer.get_uniform_location(self.program, "objectColor")
        self._specular_strength_loc = viewer.get_uniform_location(self.program, "specularStrength")
        self._ambient_strength_loc = viewer.get_uniform_location(self.program, "ambientStrength")

        self.light_position = (-20.0, 0.0, 0.0)
        self.light_color = (1.0, 1.0, 1.0)
        self.specular_strength = 0.3
        self.ambient_strength = 0.5

    def use(self):
        """Bind this shader program."""
        self.viewer.use_program(self.program)

    def set_camera(self, camera: Camera):
        """Set projection-view matrix and camera position."""
        view = camera.view_matrix()
        proj = camera.projection_matrix(0.01, 100.0)
        pv = proj @ view

        self.viewer.uniform_matrix_4fv(self._projection_view_loc, False, _mat_to_col_major(pv))
        self.viewer.uniform_3f(self._view_pos_loc, *camera.position)
        self.viewer.uniform_3f(self._light_pos_loc, *self.light_position)
        self.viewer.uniform_3f(self._light_color_loc, *self.light_color)
        self.viewer.uniform_1f(self._specular_strength_loc, self.specular_strength)
        self.viewer.uniform_1f(self._ambient_strength_loc, self.ambient_strength)

    def set_model_matrix(self, model_matrix: np.ndarray):
        """Set model transformation matrix."""
        self.viewer.uniform_matrix_4fv(self._model_loc, False, _mat_to_col_major(model_matrix))

    def set_color(self, color: tuple[float, float, float] | list[float] | np.ndarray):
        """Set object diffuse color."""
        r, g, b = float(color[0]), float(color[1]), float(color[2])
        self.viewer.uniform_3f(self._object_color_loc, r, g, b)


class CellShaderWebGL:
    """WebGL2 shader for unit cell and supercell wireframes."""

    def __init__(self, viewer):
        self.viewer = viewer
        self.program = compile_webgl_program(viewer, SIMPLE_VERTEX_SHADER, CELL_FRAGMENT_SHADER)

        self._model_loc = viewer.get_uniform_location(self.program, "model")
        self._projection_view_loc = viewer.get_uniform_location(self.program, "projectionView")
        self._cell_color_loc = viewer.get_uniform_location(self.program, "cellColor")

    def use(self):
        """Bind this shader program."""
        self.viewer.use_program(self.program)

    def set_camera(self, camera: Camera):
        """Set projection-view matrix."""
        view = camera.view_matrix()
        proj = camera.projection_matrix(0.01, 100.0)
        pv = proj @ view
        self.viewer.uniform_matrix_4fv(self._projection_view_loc, False, _mat_to_col_major(pv))

    def set_model_matrix(self, model_matrix: np.ndarray):
        """Set model transformation matrix."""
        self.viewer.uniform_matrix_4fv(self._model_loc, False, _mat_to_col_major(model_matrix))

    def set_color(self, color: tuple[float, float, float] = (0.6, 0.6, 0.6)):
        """Set wireframe color."""
        r, g, b = float(color[0]), float(color[1]), float(color[2])
        self.viewer.uniform_3f(self._cell_color_loc, r, g, b)


class AnisotropyShaderWebGL:
    """WebGL2 shader for semi-transparent anisotropy ellipsoids."""

    def __init__(self, viewer):
        self.viewer = viewer
        self.program = compile_webgl_program(viewer, OBJECT_VERTEX_SHADER, FRESNEL_FRAGMENT_SHADER)

        self._model_loc = viewer.get_uniform_location(self.program, "model")
        self._projection_view_loc = viewer.get_uniform_location(self.program, "projectionView")
        self._view_pos_loc = viewer.get_uniform_location(self.program, "viewPos")
        self._object_color_loc = viewer.get_uniform_location(self.program, "objectColor")

    def use(self):
        """Bind this shader program."""
        self.viewer.use_program(self.program)

    def set_camera(self, camera: Camera):
        """Set projection-view matrix and view position."""
        view = camera.view_matrix()
        proj = camera.projection_matrix(0.01, 100.0)
        pv = proj @ view
        self.viewer.uniform_matrix_4fv(self._projection_view_loc, False, _mat_to_col_major(pv))
        self.viewer.uniform_3f(self._view_pos_loc, *camera.position)

    def set_model_matrix(self, model_matrix: np.ndarray):
        """Set model transformation matrix."""
        self.viewer.uniform_matrix_4fv(self._model_loc, False, _mat_to_col_major(model_matrix))

    def set_color(self, color: tuple[float, float, float] | list[float] | np.ndarray):
        """Set anisotropy color."""
        r, g, b = float(color[0]), float(color[1]), float(color[2])
        self.viewer.uniform_3f(self._object_color_loc, r, g, b)


class AxesShaderWebGL:
    """WebGL2 shader for Cartesian and lattice coordinate axes."""

    def __init__(self, viewer):
        self.viewer = viewer
        self.program = compile_webgl_program(viewer, SIMPLE_VERTEX_SHADER, EMISSION_FRAGMENT_SHADER)

        self._model_loc = viewer.get_uniform_location(self.program, "model")
        self._projection_view_loc = viewer.get_uniform_location(self.program, "projectionView")
        self._object_color_loc = viewer.get_uniform_location(self.program, "objectColor")

    def use(self):
        """Bind this shader program."""
        self.viewer.use_program(self.program)

    def set_camera(self, camera: Camera):
        """Set axes orientation view matrix."""
        view = camera.axes_view_matrix()
        proj = camera.orthographic_matrix(0.01, 100.0)
        pv = proj @ view
        self.viewer.uniform_matrix_4fv(self._projection_view_loc, False, _mat_to_col_major(pv))

    def set_model_matrix(self, model_matrix: np.ndarray):
        """Set model transformation matrix."""
        self.viewer.uniform_matrix_4fv(self._model_loc, False, _mat_to_col_major(model_matrix))

    def set_color(self, color: tuple[float, float, float] | list[float] | np.ndarray):
        """Set axes arrow color."""
        r, g, b = float(color[0]), float(color[1]), float(color[2])
        self.viewer.uniform_3f(self._object_color_loc, r, g, b)


class SelectionShaderWebGL:
    """WebGL2 shader for selection outlines."""

    def __init__(self, viewer):
        self.viewer = viewer
        self.program = compile_webgl_program(viewer, SIMPLE_VERTEX_SHADER, EMISSION_FRAGMENT_SHADER)

        self._model_loc = viewer.get_uniform_location(self.program, "model")
        self._projection_view_loc = viewer.get_uniform_location(self.program, "projectionView")
        self._object_color_loc = viewer.get_uniform_location(self.program, "objectColor")

    def use(self):
        """Bind this shader program."""
        self.viewer.use_program(self.program)

    def set_camera(self, camera: Camera):
        """Set projection-view matrix."""
        view = camera.view_matrix()
        proj = camera.projection_matrix(0.01, 100.0)
        pv = proj @ view
        self.viewer.uniform_matrix_4fv(self._projection_view_loc, False, _mat_to_col_major(pv))

    def set_model_matrix(self, model_matrix: np.ndarray):
        """Set model transformation matrix."""
        self.viewer.uniform_matrix_4fv(self._model_loc, False, _mat_to_col_major(model_matrix))

    def set_color(self, color: tuple[float, float, float]):
        """Set highlight color."""
        r, g, b = float(color[0]), float(color[1]), float(color[2])
        self.viewer.uniform_3f(self._object_color_loc, r, g, b)

