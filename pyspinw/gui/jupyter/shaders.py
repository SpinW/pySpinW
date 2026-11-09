"""WebGL 2.0 (GLSL ES 3.00) shaders and ipywebgl uniform helpers for the Jupyter viewer."""

import numpy as np

from pyspinw.gui.rendering.selectionmode import SelectionMode

# ---------------------------------------------------------------------------
# GLSL ES 3.00 Vertex Shaders
# ---------------------------------------------------------------------------

OBJECT_VERTEX_SHADER = """#version 300 es
precision highp float;

layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;

layout (std140) uniform ViewBlock {
    mat4 u_cameraMatrix;
    mat4 u_viewMatrix;
    mat4 u_projectionMatrix;
    mat4 u_viewProjectionMatrix;
};

uniform mat4 model;
uniform int useCustomProjection;
uniform mat4 customProjectionMatrix;

out vec3 FragPos;
out vec3 Normal;
out vec3 ViewPos;

void main()
{
    FragPos = vec3(model * vec4(aPos, 1.0));
    Normal = mat3(transpose(inverse(model))) * aNormal;
    ViewPos = u_cameraMatrix[3].xyz;

    mat4 projView = (useCustomProjection != 0)
        ? (customProjectionMatrix * u_viewMatrix)
        : u_viewProjectionMatrix;

    gl_Position = projView * vec4(FragPos, 1.0);
}
"""

SIMPLE_VERTEX_SHADER = """#version 300 es
precision highp float;

layout (location = 0) in vec3 position;
layout (location = 1) in vec3 aNormal;

layout (std140) uniform ViewBlock {
    mat4 u_cameraMatrix;
    mat4 u_viewMatrix;
    mat4 u_projectionMatrix;
    mat4 u_viewProjectionMatrix;
};

uniform mat4 model;
uniform int useCustomProjection;
uniform mat4 customProjectionMatrix;

void main()
{
    mat4 projView = (useCustomProjection != 0)
        ? (customProjectionMatrix * u_viewMatrix)
        : u_viewProjectionMatrix;

    gl_Position = projView * model * vec4(position, 1.0);
}
"""

AXES_VERTEX_SHADER = """#version 300 es
precision highp float;

layout (location = 0) in vec3 position;
layout (location = 1) in vec3 aNormal;

layout (std140) uniform ViewBlock {
    mat4 u_cameraMatrix;
    mat4 u_viewMatrix;
    mat4 u_projectionMatrix;
    mat4 u_viewProjectionMatrix;
};

uniform mat4 model;

void main()
{
    vec3 worldPos = vec3(model * vec4(position, 1.0));
    vec3 viewPos = mat3(u_viewMatrix) * worldPos;
    gl_Position = vec4(viewPos.x / 1.5, viewPos.y / 1.5, -viewPos.z / 10.0, 1.0);
}
"""

# ---------------------------------------------------------------------------
# GLSL ES 3.00 Fragment Shaders
# ---------------------------------------------------------------------------

OBJECT_FRAGMENT_SHADER = """#version 300 es
precision highp float;

in vec3 FragPos;
in vec3 Normal;
in vec3 ViewPos;

out vec4 FragColor;

uniform vec3 lightPos;
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

    vec3 viewDir = normalize(ViewPos - FragPos);
    vec3 halfwayDir = normalize(lightDir + viewDir);
    float spec = pow(max(dot(norm, halfwayDir), 0.0), 32.0);

    vec3 specular = specularStrength * spec * lightColor;

    vec3 result = (ambient + diffuse) * objectColor + specular;

    FragColor = vec4(result, 1.0);
}
"""

EMISSION_FRAGMENT_SHADER = """#version 300 es
precision highp float;

out vec4 FragColor;

uniform vec3 objectColor;

void main()
{
    FragColor = vec4(objectColor, 1.0);
}
"""

FRESNEL_FRAGMENT_SHADER = """#version 300 es
precision highp float;

in vec3 FragPos;
in vec3 Normal;
in vec3 ViewPos;

out vec4 FragColor;

uniform vec3 objectColor;

void main()
{
    vec3 N = normalize(Normal);
    vec3 V = normalize(ViewPos - FragPos);

    float alpha = pow(1.0 - max(dot(N, V), 0.0), 2.0);

    FragColor = vec4(objectColor, alpha);
}
"""

CELL_FRAGMENT_SHADER = """#version 300 es
precision highp float;

out vec4 FragColor;

void main()
{
    FragColor = vec4(0.6, 0.6, 0.6, 1.0);
}
"""


# ---------------------------------------------------------------------------
# Buffer / Uniform Formatting Helpers for ipywebgl
# ---------------------------------------------------------------------------

def to_webgl_matrix(matrix: np.ndarray) -> np.ndarray:
    """Convert a row-major (4, 4) NumPy matrix into a column-major float32 array for WebGL2.

    WebGL 2 requires ``transpose=False`` in ``gl.uniformMatrix4fv``, so row-major
    NumPy transformation matrices must be transposed before serialization.

    Parameters
    ----------
    matrix : (4, 4) array_like
        Row-major 4x4 transformation matrix.

    Returns
    -------
    (4, 4) numpy.ndarray
        Contiguous ``float32`` transposed matrix ready for ``GLViewer.uniform_matrix``.
    """
    return np.ascontiguousarray(np.asarray(matrix, dtype=np.float32).T)


def to_webgl_vec3(vec) -> np.ndarray:
    """Convert a 3-element sequence into a contiguous 1D float32 array for ``GLViewer.uniform``.

    Parameters
    ----------
    vec : array_like
        3-element vector or RGB tuple.

    Returns
    -------
    (3,) numpy.ndarray
        1D ``float32`` array of length 3.
    """
    return np.ascontiguousarray(np.asarray(vec, dtype=np.float32).reshape(3))


def to_webgl_float(value: float) -> np.ndarray:
    """Convert a scalar float into a 1-element float32 array for ``GLViewer.uniform``.

    Parameters
    ----------
    value : float
        Scalar float value.

    Returns
    -------
    (1,) numpy.ndarray
        1D ``float32`` array of length 1.
    """
    return np.array([value], dtype=np.float32)


def to_webgl_int(value: int) -> np.ndarray:
    """Convert a scalar int into a 1-element int32 array for ``GLViewer.uniform``.

    Parameters
    ----------
    value : int
        Scalar integer value.

    Returns
    -------
    (1,) numpy.ndarray
        1D ``int32`` array of length 1.
    """
    return np.array([value], dtype=np.int32)


# ---------------------------------------------------------------------------
# Shader Program Wrappers for ipywebgl.GLViewer
# ---------------------------------------------------------------------------

class WebGLObjectShader:
    """Wrapper for the Blinn-Phong object shader program in ``ipywebgl``.

    Parameters
    ----------
    gl : ipywebgl.GLViewer
        Target WebGL viewer widget.
    auto_execute : bool, optional
        Whether to immediately execute program creation commands. Default is ``False``.
    """

    def __init__(self, gl, auto_execute: bool = False):
        self.gl = gl
        self.program = gl.create_program_ext(
            OBJECT_VERTEX_SHADER,
            OBJECT_FRAGMENT_SHADER,
            attribute_location={"aPos": 0, "aNormal": 1},
            auto_execute=auto_execute,
        )
        self.light_position = (-20.0, 0.0, 0.0)
        self.light_color = (1.0, 1.0, 1.0)
        self.specular_strength = 0.3
        self.ambient_strength = 0.5

    def use(
        self,
        use_custom_projection: bool = False,
        custom_projection_matrix: np.ndarray | None = None,
    ):
        """Bind the program and set scene-level lighting and projection uniforms.

        Parameters
        ----------
        use_custom_projection : bool, optional
            Whether to override ``ViewBlock`` projection with ``custom_projection_matrix``.
            Default is ``False``.
        custom_projection_matrix : (4, 4) numpy.ndarray or None, optional
            Custom row-major 4x4 projection matrix (e.g., orthographic). Default is ``None``.
        """
        self.gl.use_program(self.program)
        self.gl.uniform("lightPos", to_webgl_vec3(self.light_position))
        self.gl.uniform("lightColor", to_webgl_vec3(self.light_color))
        self.gl.uniform("specularStrength", to_webgl_float(self.specular_strength))
        self.gl.uniform("ambientStrength", to_webgl_float(self.ambient_strength))
        self.gl.uniform("useCustomProjection", to_webgl_int(1 if use_custom_projection else 0))
        if use_custom_projection and custom_projection_matrix is not None:
            self.gl.uniform_matrix("customProjectionMatrix", to_webgl_matrix(custom_projection_matrix))

    def set_object(self, model_matrix: np.ndarray, object_color: tuple[float, float, float]):
        """Set per-object model matrix and RGB color uniforms.

        Parameters
        ----------
        model_matrix : (4, 4) numpy.ndarray
            Row-major 4x4 model transformation matrix.
        object_color : tuple of (float, float, float)
            RGB color in ``[0, 1]``.
        """
        self.gl.uniform_matrix("model", to_webgl_matrix(model_matrix))
        self.gl.uniform("objectColor", to_webgl_vec3(object_color))


class WebGLSelectionShader:
    """Wrapper for the back-face selection/hover outline shader program in ``ipywebgl``.

    Parameters
    ----------
    gl : ipywebgl.GLViewer
        Target WebGL viewer widget.
    auto_execute : bool, optional
        Whether to immediately execute program creation commands. Default is ``False``.
    """

    def __init__(self, gl, auto_execute: bool = False):
        self.gl = gl
        self.program = gl.create_program_ext(
            SIMPLE_VERTEX_SHADER,
            EMISSION_FRAGMENT_SHADER,
            attribute_location={"position": 0, "aNormal": 1},
            auto_execute=auto_execute,
        )
        self.selected_color = (1.0, 0.6, 0.1)
        self.hover_color = (1.0, 1.0, 1.0)
        self.selected_hover_color = (1.0, 0.8, 0.1)

    def color_for_mode(self, mode: SelectionMode) -> tuple[float, float, float]:
        """Resolve the highlight RGB color for a given ``SelectionMode``.

        Parameters
        ----------
        mode : SelectionMode
            Selection state of the object.

        Returns
        -------
        tuple of (float, float, float)
            RGB color tuple.
        """
        match mode:
            case SelectionMode.SELECTED:
                return self.selected_color
            case SelectionMode.HOVER:
                return self.hover_color
            case SelectionMode.SELECTED_HOVER:
                return self.selected_hover_color
            case _:
                return (1.0, 1.0, 1.0)

    def use(
        self,
        use_custom_projection: bool = False,
        custom_projection_matrix: np.ndarray | None = None,
    ):
        """Bind the selection shader program and configure projection uniforms.

        Parameters
        ----------
        use_custom_projection : bool, optional
            Whether to override ``ViewBlock`` projection with ``custom_projection_matrix``.
            Default is ``False``.
        custom_projection_matrix : (4, 4) numpy.ndarray or None, optional
            Custom row-major 4x4 projection matrix. Default is ``None``.
        """
        self.gl.use_program(self.program)
        self.gl.uniform("useCustomProjection", to_webgl_int(1 if use_custom_projection else 0))
        if use_custom_projection and custom_projection_matrix is not None:
            self.gl.uniform_matrix("customProjectionMatrix", to_webgl_matrix(custom_projection_matrix))

    def set_object(self, model_matrix: np.ndarray, mode: SelectionMode):
        """Set per-object model matrix and selection mode color.

        Parameters
        ----------
        model_matrix : (4, 4) numpy.ndarray
            Row-major 4x4 model transformation matrix.
        mode : SelectionMode
            Selection state of the object.
        """
        self.gl.uniform_matrix("model", to_webgl_matrix(model_matrix))
        self.gl.uniform("objectColor", to_webgl_vec3(self.color_for_mode(mode)))


class WebGLAnisotropyShader:
    """Wrapper for the Fresnel alpha-blended anisotropy shader program in ``ipywebgl``.

    Parameters
    ----------
    gl : ipywebgl.GLViewer
        Target WebGL viewer widget.
    auto_execute : bool, optional
        Whether to immediately execute program creation commands. Default is ``False``.
    """

    def __init__(self, gl, auto_execute: bool = False):
        self.gl = gl
        self.program = gl.create_program_ext(
            OBJECT_VERTEX_SHADER,
            FRESNEL_FRAGMENT_SHADER,
            attribute_location={"aPos": 0, "aNormal": 1},
            auto_execute=auto_execute,
        )

    def use(
        self,
        object_color: tuple[float, float, float],
        use_custom_projection: bool = False,
        custom_projection_matrix: np.ndarray | None = None,
    ):
        """Bind the anisotropy shader program and set color and projection uniforms.

        Parameters
        ----------
        object_color : tuple of (float, float, float)
            RGB color for anisotropy ellipsoids.
        use_custom_projection : bool, optional
            Whether to override ``ViewBlock`` projection with ``custom_projection_matrix``.
            Default is ``False``.
        custom_projection_matrix : (4, 4) numpy.ndarray or None, optional
            Custom row-major 4x4 projection matrix. Default is ``None``.
        """
        self.gl.use_program(self.program)
        self.gl.uniform("objectColor", to_webgl_vec3(object_color))
        self.gl.uniform("useCustomProjection", to_webgl_int(1 if use_custom_projection else 0))
        if use_custom_projection and custom_projection_matrix is not None:
            self.gl.uniform_matrix("customProjectionMatrix", to_webgl_matrix(custom_projection_matrix))

    def set_model_matrix(self, model_matrix: np.ndarray):
        """Set the per-ellipsoid model matrix uniform.

        Parameters
        ----------
        model_matrix : (4, 4) numpy.ndarray
            Row-major 4x4 model transformation matrix.
        """
        self.gl.uniform_matrix("model", to_webgl_matrix(model_matrix))


class WebGLCellShader:
    """Wrapper for the unit cell / supercell wireframe shader program in ``ipywebgl``.

    Parameters
    ----------
    gl : ipywebgl.GLViewer
        Target WebGL viewer widget.
    auto_execute : bool, optional
        Whether to immediately execute program creation commands. Default is ``False``.
    """

    def __init__(self, gl, auto_execute: bool = False):
        self.gl = gl
        self.program = gl.create_program_ext(
            SIMPLE_VERTEX_SHADER,
            CELL_FRAGMENT_SHADER,
            attribute_location={"position": 0},
            auto_execute=auto_execute,
        )

    def use(
        self,
        use_custom_projection: bool = False,
        custom_projection_matrix: np.ndarray | None = None,
    ):
        """Bind the cell wireframe shader program and set projection uniforms.

        Parameters
        ----------
        use_custom_projection : bool, optional
            Whether to override ``ViewBlock`` projection with ``custom_projection_matrix``.
            Default is ``False``.
        custom_projection_matrix : (4, 4) numpy.ndarray or None, optional
            Custom row-major 4x4 projection matrix. Default is ``None``.
        """
        self.gl.use_program(self.program)
        self.gl.uniform("useCustomProjection", to_webgl_int(1 if use_custom_projection else 0))
        if use_custom_projection and custom_projection_matrix is not None:
            self.gl.uniform_matrix("customProjectionMatrix", to_webgl_matrix(custom_projection_matrix))

    def set_model_matrix(self, model_matrix: np.ndarray):
        """Set the wireframe cube model matrix uniform.

        Parameters
        ----------
        model_matrix : (4, 4) numpy.ndarray
            Row-major 4x4 model transformation matrix.
        """
        self.gl.uniform_matrix("model", to_webgl_matrix(model_matrix))


class WebGLAxesShader:
    """Wrapper for the corner orientation axes shader program in ``ipywebgl``.

    Parameters
    ----------
    gl : ipywebgl.GLViewer
        Target WebGL viewer widget.
    auto_execute : bool, optional
        Whether to immediately execute program creation commands. Default is ``False``.
    """

    def __init__(self, gl, auto_execute: bool = False):
        self.gl = gl
        self.program = gl.create_program_ext(
            AXES_VERTEX_SHADER,
            EMISSION_FRAGMENT_SHADER,
            attribute_location={"position": 0, "aNormal": 1},
            auto_execute=auto_execute,
        )

    def use(self):
        """Bind the axes shader program."""
        self.gl.use_program(self.program)

    def set_axis(self, model_matrix: np.ndarray, object_color: tuple[float, float, float]):
        """Set the model matrix and color for a single axis arrow.

        Parameters
        ----------
        model_matrix : (4, 4) numpy.ndarray
            Row-major 4x4 model transformation matrix.
        object_color : tuple of (float, float, float)
            RGB color tuple for the axis arrow.
        """
        self.gl.uniform_matrix("model", to_webgl_matrix(model_matrix))
        self.gl.uniform("objectColor", to_webgl_vec3(object_color))
