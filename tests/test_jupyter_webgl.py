"""Unit tests for Jupyter WebGL shaders and model geometry generators."""

import numpy as np
import pytest

from pyspinw.gui.jupyter import check_jupyter_dependencies
from pyspinw.gui.jupyter.models_webgl import (
    generate_arrow_geometry,
    generate_cube_wireframe_geometry,
    generate_sphere_geometry,
    generate_tube_geometry,
    WebGLMesh,
)
from pyspinw.gui.jupyter.shaders_webgl import (
    AnisotropyShaderWebGL,
    AxesShaderWebGL,
    CELL_FRAGMENT_SHADER,
    CellShaderWebGL,
    EMISSION_FRAGMENT_SHADER,
    FRESNEL_FRAGMENT_SHADER,
    OBJECT_FRAGMENT_SHADER,
    OBJECT_VERTEX_SHADER,
    ObjectShaderWebGL,
    SIMPLE_VERTEX_SHADER,
    SelectionShaderWebGL,
    _mat_to_col_major,
)
from pyspinw.gui.camera import Camera


class MockGLViewer:
    """Mock ipywebgl.GLViewer for unit testing shaders and mesh binding."""

    def __init__(self):
        self._next_id = 1
        self.calls = []

    def create_shader(self, shader_type):
        sid = f"shader_{self._next_id}"
        self._next_id += 1
        self.calls.append(("create_shader", shader_type, sid))
        return sid

    def shader_source(self, shader, source):
        self.calls.append(("shader_source", shader, len(source)))

    def compile_shader(self, shader):
        self.calls.append(("compile_shader", shader))

    def create_program(self):
        pid = f"prog_{self._next_id}"
        self._next_id += 1
        self.calls.append(("create_program", pid))
        return pid

    def attach_shader(self, prog, shader):
        self.calls.append(("attach_shader", prog, shader))

    def link_program(self, prog):
        self.calls.append(("link_program", prog))

    def use_program(self, prog):
        self.calls.append(("use_program", prog))

    def uniform_matrix(self, name, array):
        self.calls.append(("uniform_matrix", name, len(array)))

    def uniform(self, name, array):
        self.calls.append(("uniform", name, len(array)))

    def get_uniform_location(self, prog, name):
        return f"loc_{name}"

    def uniform_matrix_4fv(self, loc, transpose, data):
        self.calls.append(("uniform_matrix_4fv", loc, transpose, len(data)))

    def uniform_3f(self, loc, x, y, z):
        self.calls.append(("uniform_3f", loc, x, y, z))

    def uniform_1f(self, loc, val):
        self.calls.append(("uniform_1f", loc, val))

    def create_buffer(self):
        bid = f"buf_{self._next_id}"
        self._next_id += 1
        return bid

    def bind_buffer(self, target, buf):
        self.calls.append(("bind_buffer", target, buf))

    def buffer_data(self, target, data, usage):
        self.calls.append(("buffer_data", target, len(data), usage))

    def create_vertex_array(self):
        aid = f"vao_{self._next_id}"
        self._next_id += 1
        return aid

    def bind_vertex_array(self, vao):
        self.calls.append(("bind_vertex_array", vao))

    def enable_vertex_attrib_array(self, index):
        self.calls.append(("enable_vertex_attrib_array", index))

    def vertex_attrib_pointer(self, index, size, type_, norm, stride, offset):
        self.calls.append(("vertex_attrib_pointer", index, size, type_, stride, offset))

    def enable(self, cap):
        self.calls.append(("enable", cap))

    def cull_face(self, mode):
        self.calls.append(("cull_face", mode))

    def draw_arrays(self, mode, first, count):
        self.calls.append(("draw_arrays", mode, first, count))


def test_dependency_checker():
    """Verify check_jupyter_dependencies raises ImportError if packages missing."""
    # When tested in an environment where ipywebgl/ipyevents might not be installed,
    # it should raise an informative ImportError.
    try:
        check_jupyter_dependencies()
    except ImportError as e:
        assert "pyspinw[jupyter]" in str(e)


def test_glsl_sources():
    """Ensure GLSL ES 3.00 shader sources are valid version strings and non-empty."""
    for src in [OBJECT_VERTEX_SHADER, OBJECT_FRAGMENT_SHADER, SIMPLE_VERTEX_SHADER,
                CELL_FRAGMENT_SHADER, FRESNEL_FRAGMENT_SHADER, EMISSION_FRAGMENT_SHADER]:
        assert src.strip().startswith("#version 300 es")
        assert len(src) > 50


def test_mat_to_col_major():
    """Verify that matrix conversion produces flat 16-element float32 column-major arrays."""
    mat = np.array([
        [1, 2, 3, 4],
        [5, 6, 7, 8],
        [9, 10, 11, 12],
        [13, 14, 15, 16],
    ], dtype=np.float32)

    flat = _mat_to_col_major(mat)
    assert len(flat) == 16
    assert flat.dtype == np.float32
    # In column-major, first 4 entries are column 0: 1, 5, 9, 13
    np.testing.assert_array_equal(flat[:4], [1, 5, 9, 13])


def test_sphere_geometry_generation():
    """Verify sphere vertex and normal generation."""
    sphere = generate_sphere_geometry(divisions=2, radius=1.0)
    assert sphere.ndim == 2
    assert sphere.shape[1] == 6
    assert sphere.dtype == np.float32
    assert np.all(np.isfinite(sphere))

    # Test with padding
    sphere_pad = generate_sphere_geometry(divisions=2, radius=0.1, padding=0.02)
    assert sphere_pad.shape[1] == 6
    assert len(sphere_pad) == len(sphere)


def test_tube_geometry_generation():
    """Verify tube vertex and normal generation."""
    tube = generate_tube_geometry(n_points=10)
    assert tube.ndim == 2
    assert tube.shape[1] == 6
    assert tube.dtype == np.float32
    assert len(tube) == (10 - 1) * 6
    assert np.all(np.isfinite(tube))


def test_arrow_geometry_generation():
    """Verify arrow vertex and normal generation."""
    arrow = generate_arrow_geometry(n_points=10)
    assert arrow.ndim == 2
    assert arrow.shape[1] == 6
    assert arrow.dtype == np.float32
    assert np.all(np.isfinite(arrow))

    # Centered at z in range [-0.5, 0.5]
    z_coords = arrow[:, 2]
    assert np.isclose(z_coords.min(), -0.5, atol=0.1)
    assert np.isclose(z_coords.max(), 0.5, atol=0.1)


def test_cube_wireframe_generation():
    """Verify cube wireframe vertex generation."""
    cube = generate_cube_wireframe_geometry()
    assert cube.shape == (24, 3)
    assert cube.dtype == np.float32
    assert np.all(np.isfinite(cube))
    assert cube.min() == 0.0
    assert cube.max() == 1.0


def test_webgl_shaders_with_mock():
    """Verify that WebGL shader wrappers properly initialize and bind uniforms."""
    viewer = MockGLViewer()
    cam = Camera()

    # Object shader
    obj_sh = ObjectShaderWebGL(viewer)
    obj_sh.use()
    obj_sh.set_camera(cam)
    obj_sh.set_model_matrix(np.eye(4))
    obj_sh.set_color((1.0, 0.0, 0.0))

    # Cell shader
    cell_sh = CellShaderWebGL(viewer)
    cell_sh.use()
    cell_sh.set_camera(cam)
    cell_sh.set_model_matrix(np.eye(4))
    cell_sh.set_color((0.5, 0.5, 0.5))

    # Anisotropy shader
    ani_sh = AnisotropyShaderWebGL(viewer)
    ani_sh.use()
    ani_sh.set_camera(cam)
    ani_sh.set_model_matrix(np.eye(4))
    ani_sh.set_color((0.2, 0.5, 1.0))

    # Axes shader
    axes_sh = AxesShaderWebGL(viewer)
    axes_sh.use()
    axes_sh.set_camera(cam)
    axes_sh.set_model_matrix(np.eye(4))
    axes_sh.set_color((0.8, 0.4, 0.0))

    # Selection shader
    sel_sh = SelectionShaderWebGL(viewer)
    sel_sh.use()
    sel_sh.set_camera(cam)
    sel_sh.set_model_matrix(np.eye(4))
    sel_sh.set_color((1.0, 0.8, 0.1))

    assert len(viewer.calls) > 20


def test_webgl_mesh_with_mock():
    """Verify that WebGLMesh properly registers VBOs and VAOs and issues draw calls."""
    viewer = MockGLViewer()
    verts = generate_sphere_geometry(divisions=2)

    mesh = WebGLMesh(viewer, verts, has_normals=True, draw_mode="TRIANGLES")
    assert mesh.vertex_count == len(verts)

    mesh.render_triangles()
    mesh.render_back()

    wireframe = generate_cube_wireframe_geometry()
    wireframe_mesh = WebGLMesh(viewer, wireframe, has_normals=False, draw_mode="LINES")
    assert wireframe_mesh.vertex_count == 24
    wireframe_mesh.render_wireframe()

    # Check draw_arrays calls
    draw_calls = [c for c in viewer.calls if c[0] == "draw_arrays"]
    assert len(draw_calls) == 3

