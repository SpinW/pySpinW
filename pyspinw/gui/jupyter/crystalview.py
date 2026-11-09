"""Jupyter WebGL crystal structure viewer widget built on ``ipywebgl``."""

from collections.abc import Callable

import ipywebgl
import numpy as np

from pyspinw.gui.camera import Camera
from pyspinw.gui.displayoptions import DisplayOptions
from pyspinw.gui.jupyter.shaders import (
    WebGLAnisotropyShader,
    WebGLAxesShader,
    WebGLCellShader,
    WebGLObjectShader,
    WebGLSelectionShader,
)
from pyspinw.gui.rendermodel import RenderModel
from pyspinw.gui.rendering.models.arrow import Arrow
from pyspinw.gui.rendering.models.sphere import Sphere
from pyspinw.gui.rendering.models.tube import Tube
from pyspinw.gui.rendering.models.wrireframe_cube import WireframeCube
from pyspinw.gui.rendering.selectionmode import SelectionMode

# Transforms for xyz arrows in axes viewport
_z_axis_mat = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0.5], [0, 0, 0, 1]], dtype=np.float32)
_y_axis_mat = np.array([[1, 0, 0, 0], [0, 0, 1, 0.5], [0, -1, 0, 0], [0, 0, 0, 1]], dtype=np.float32)
_x_axis_mat = np.array([[0, 0, 1, 0.5], [0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1]], dtype=np.float32)

_selection_size = 0.2


class WebGLMesh:
    """WebGL vertex buffer and vertex array object wrapper for triangle meshes.

    Parameters
    ----------
    gl : ipywebgl.GLViewer
        Target WebGL viewer widget.
    program : ipywebgl.glresource.GLResourceWidget
        Reference shader program used when binding vertex attributes.
    vertices_and_normals : (n, 6) numpy.ndarray or list of (n, 6) numpy.ndarray
        Interleaved ``float32`` 3D vertex positions and 3D normals.
    """

    def __init__(self, gl: ipywebgl.GLViewer, program, vertices_and_normals: np.ndarray | list[np.ndarray]):
        if isinstance(vertices_and_normals, list):
            data = np.ascontiguousarray(np.vstack(vertices_and_normals), dtype=np.float32)
        else:
            data = np.ascontiguousarray(vertices_and_normals.reshape(-1, 6), dtype=np.float32)

        self.vertex_count = int(data.shape[0])
        self.vbo = gl.create_buffer_ext(
            target="ARRAY_BUFFER",
            src_data=data,
            usage="STATIC_DRAW",
            auto_execute=False,
        )
        self.vao = gl.create_vertex_array_ext(
            program,
            [(self.vbo, "3f32 3f32", 0, 1)],
            auto_execute=False,
        )

    def render_triangles(self, gl: ipywebgl.GLViewer):
        """Append front-face triangle draw commands to ``gl``.

        Parameters
        ----------
        gl : ipywebgl.GLViewer
            Target WebGL viewer widget.
        """
        gl.enable(cull_face=True)
        gl.cull_face("BACK")
        gl.bind_vertex_array(self.vao)
        gl.draw_arrays("TRIANGLES", 0, self.vertex_count)
        gl.bind_vertex_array(None)

    def render_back(self, gl: ipywebgl.GLViewer):
        """Append back-face triangle draw commands to ``gl`` for selection outlines.

        Parameters
        ----------
        gl : ipywebgl.GLViewer
            Target WebGL viewer widget.
        """
        gl.enable(cull_face=True)
        gl.cull_face("FRONT")
        gl.bind_vertex_array(self.vao)
        gl.draw_arrays("TRIANGLES", 0, self.vertex_count)
        gl.bind_vertex_array(None)


class WebGLWireframeMesh:
    """WebGL vertex buffer and vertex array object wrapper for wireframe line meshes.

    Parameters
    ----------
    gl : ipywebgl.GLViewer
        Target WebGL viewer widget.
    program : ipywebgl.glresource.GLResourceWidget
        Reference shader program used when binding vertex attributes.
    vertices : (n, 3) numpy.ndarray
        ``float32`` array of 3D line endpoint coordinates.
    """

    def __init__(self, gl: ipywebgl.GLViewer, program, vertices: np.ndarray):
        data = np.ascontiguousarray(vertices.reshape(-1, 3), dtype=np.float32)
        self.vertex_count = int(data.shape[0])
        self.vbo = gl.create_buffer_ext(
            target="ARRAY_BUFFER",
            src_data=data,
            usage="STATIC_DRAW",
            auto_execute=False,
        )
        self.vao = gl.create_vertex_array_ext(
            program,
            [(self.vbo, "3f32", 0)],
            auto_execute=False,
        )

    def render_wireframe(self, gl: ipywebgl.GLViewer):
        """Append ``LINES`` draw commands to ``gl``.

        Parameters
        ----------
        gl : ipywebgl.GLViewer
            Target WebGL viewer widget.
        """
        gl.bind_vertex_array(self.vao)
        gl.draw_arrays("LINES", 0, self.vertex_count)
        gl.bind_vertex_array(None)


class JupyterCrystalView:
    """WebGL 3D crystal structure and magnetic interaction view for Jupyter notebooks.

    Parameters
    ----------
    render_model : RenderModel
        Model containing the crystal structure and magnetic interactions to render.
    initial_rotation : (3, 3) numpy.ndarray or None, optional
        Initial 3x3 view rotation matrix. Default is ``None`` (identity matrix).
    initial_distance : float or None, optional
        Initial camera distance from the view origin. Default is ``None`` (``10.0``).
    width : int, optional
        Canvas width in pixels. Default is ``640``.
    height : int, optional
        Canvas height in pixels. Default is ``520``.
    """

    axes_size = 100
    axes_padding = 10

    def __init__(
        self,
        render_model: RenderModel,
        initial_rotation: np.ndarray | None = None,
        initial_distance: float | None = None,
        width: int = 640,
        height: int = 520,
    ):
        self.render_model = render_model
        self._display_options = DisplayOptions()
        self.camera = Camera(horizontal_pixels=width, vertical_pixels=height)

        self._initial_rotation = np.eye(3, dtype=float) if initial_rotation is None else np.array(initial_rotation)
        self._initial_distance = 10.0 if initial_distance is None else float(initial_distance)

        self.view_origin = np.array(render_model.expanded.structure.unit_cell.centre, dtype=float)
        self.view_rotation = self._initial_rotation.copy()
        self.view_radius = self._initial_distance

        self._hover_ids: list[int] = []
        self._current_selection: list[int] = []

        self._hover_callbacks: list[Callable[[], None]] = []
        self._selection_callbacks: list[Callable[[], None]] = []

        self.gl = ipywebgl.GLViewer(
            width=width,
            height=height,
            camera_fov=50.0,
            camera_near=0.05,
            camera_far=500.0,
            mouse_speed=1.0,
            move_speed=0.25,
            shader_matrix_major="row_major",
        )

        self._apply_camera_pose(self.view_origin, self.view_rotation, self.view_radius)
        self._init_gl_resources()
        self.update_scene()

    @property
    def widget(self) -> ipywebgl.GLViewer:
        """Underlying ``ipywebgl.GLViewer`` widget."""
        return self.gl

    @property
    def display_options(self) -> DisplayOptions:
        """Current display options."""
        return self._display_options

    @display_options.setter
    def display_options(self, value: DisplayOptions):
        self._display_options = value
        self.update_scene()

    @property
    def hover_ids(self) -> list[int]:
        """Render IDs of currently hovered elements."""
        return self._hover_ids

    @hover_ids.setter
    def hover_ids(self, value: list[int]):
        self._hover_ids = list(value)
        self.update_scene()

    @property
    def current_selection(self) -> list[int]:
        """Render IDs of currently selected elements."""
        return self._current_selection

    @current_selection.setter
    def current_selection(self, value: list[int]):
        self._current_selection = list(value)
        self.update_scene()

    def on_hover_changed(self, callback: Callable[[], None]):
        """Register a callback invoked when hovered render IDs change.

        Parameters
        ----------
        callback : callable
            Zero-argument function to call on hover changes.
        """
        self._hover_callbacks.append(callback)

    def on_selection_changed(self, callback: Callable[[], None]):
        """Register a callback invoked when selected render IDs change.

        Parameters
        ----------
        callback : callable
            Zero-argument function to call on selection changes.
        """
        self._selection_callbacks.append(callback)

    def _apply_camera_pose(self, origin: np.ndarray, rotation: np.ndarray, radius: float):
        """Compute and assign ``camera_pos``, ``camera_yaw``, and ``camera_pitch`` on ``self.gl``.

        Parameters
        ----------
        origin : (3,) numpy.ndarray
            World-space target origin.
        rotation : (3, 3) numpy.ndarray
            3x3 camera orientation matrix.
        radius : float
            Distance from ``origin`` to the camera.
        """
        forward_from_target = rotation @ np.array([0.0, 0.0, 1.0], dtype=float)
        norm = np.linalg.norm(forward_from_target)
        if norm > 1e-12:
            forward_from_target = forward_from_target / norm

        pitch_rad = -np.arcsin(np.clip(forward_from_target[1], -1.0, 1.0))
        yaw_rad = np.arctan2(forward_from_target[0], forward_from_target[2])

        camera_world = origin + rotation @ np.array([0.0, 0.0, radius], dtype=float)

        self.gl.camera_pos = [float(camera_world[0]), float(camera_world[1]), float(camera_world[2])]
        self.gl.camera_yaw = float(np.degrees(yaw_rad))
        self.gl.camera_pitch = float(np.degrees(pitch_rad))

        self.camera.position = tuple(float(x) for x in camera_world)
        self.camera.look_at = tuple(float(x) for x in origin)
        self.camera.up = tuple(float(x) for x in (rotation @ np.array([0.0, 1.0, 0.0], dtype=float)))

    def _init_gl_resources(self):
        """Create and upload shaders, VBOs, and VAOs to the WebGL context."""
        self.object_shader = WebGLObjectShader(self.gl, auto_execute=False)
        self.selection_shader = WebGLSelectionShader(self.gl, auto_execute=False)
        self.anisotropy_shader = WebGLAnisotropyShader(self.gl, auto_execute=False)
        self.cell_shader = WebGLCellShader(self.gl, auto_execute=False)
        self.axes_shader = WebGLAxesShader(self.gl, auto_execute=False)

        ref_prog = self.object_shader.program

        self.sphere = WebGLMesh(self.gl, ref_prog, Sphere.build_geometry(3))
        self.selection_sphere = WebGLMesh(self.gl, ref_prog, Sphere.build_geometry(3, padding=0.1))

        self.small_sphere = WebGLMesh(self.gl, ref_prog, Sphere.build_geometry(3, radius=0.1))
        self.selection_small_sphere = WebGLMesh(self.gl, ref_prog, Sphere.build_geometry(3, radius=0.1, padding=0.01))
        self.anisotropy_sphere = WebGLMesh(self.gl, ref_prog, Sphere.build_geometry(5, radius=0.1))

        self.tube = WebGLMesh(self.gl, ref_prog, Tube.build_geometry())
        self.selection_tube = WebGLMesh(self.gl, ref_prog, Tube.build_geometry(padding=0.2))

        self.arrow = WebGLMesh(self.gl, ref_prog, Arrow.build_geometry())
        self.selection_arrow = WebGLMesh(self.gl, ref_prog, Arrow.build_geometry(padding=0.03, modify_ring=True))

        self.cube = WebGLWireframeMesh(self.gl, self.cell_shader.program, WireframeCube.build_geometry())

        # Execute one-time resource creation commands
        self.gl.execute_commands(execute_once=True, clear_previous=True)

    def reset_view(self):
        """Reset the camera view to the default origin, rotation, and distance."""
        self.view_origin = np.array(self.render_model.expanded.structure.unit_cell.centre, dtype=float)
        self.view_rotation = np.eye(3, dtype=float)
        self.view_radius = 10.0
        self._apply_camera_pose(self.view_origin, self.view_rotation, self.view_radius)
        self.update_scene()

    def click_on_element(self, index: int, shift: bool = False):
        """Update selection state as if clicking on an element with ``render_id == index``.

        Parameters
        ----------
        index : int
            Render ID of the clicked element (``0`` if background).
        shift : bool, optional
            Whether to append to the current selection instead of replacing it.
            Default is ``False``.
        """
        if index == 0:
            if not shift:
                self._current_selection = []
                self.update_scene()
                for cb in self._selection_callbacks:
                    cb()
        else:
            if not shift:
                self._current_selection = [index]
            elif index not in self._current_selection:
                self._current_selection.append(index)
            self.update_scene()
            for cb in self._selection_callbacks:
                cb()

    def update_scene(self):
        """Build and dispatch the WebGL rendering command buffer for the current scene state."""
        w, h = int(self.gl.width), int(self.gl.height)
        self.camera.horizontal_pixels = w
        self.camera.vertical_pixels = h
        self.camera.position = tuple(float(x) for x in self.gl.camera_pos)
        self.camera.look_at = tuple(float(x) for x in self.view_origin)
        self.camera.orthographic = not self._display_options.perspective
        self.camera.fixed_size = None

        use_custom_proj = not self._display_options.perspective
        custom_proj = self.camera.orthographic_matrix(0.05, 500.0) if use_custom_proj else None

        self.gl.viewport(0, 0, w, h)
        self.gl.enable(depth_test=True)
        self.gl.enable(cull_face=True)
        self.gl.disable(blend=True)

        bg_r, bg_g, bg_b = self._display_options.background_color
        self.gl.clear_color(bg_r, bg_g, bg_b, 1.0)
        self.gl.clear(color_bit_buffer=True, depth_buffer_bit=True, stencil_buffer_bit=False)

        # Configure selection highlight colors
        self.selection_shader.selected_color = self._display_options.selected_color
        self.selection_shader.hover_color = self._display_options.hover_color
        self.selection_shader.selected_hover_color = self._display_options.selected_hover_color

        spin_scale = 2.0 * self._display_options.atom_spin_scaling
        spin_scale_matrix = np.diag([spin_scale, spin_scale, spin_scale, 1.0]).astype(np.float32)

        anisotropy_scale = 3.0 * self._display_options.anisotropy_scaling
        anisotropy_scale_matrix = np.diag([anisotropy_scale, anisotropy_scale, anisotropy_scale, 1.0]).astype(
            np.float32
        )

        exchange_scale = 0.1 * self._display_options.exchange_scaling
        exchange_scaling = np.diag([exchange_scale, exchange_scale, 1.0, 1.0]).astype(np.float32)
        exchange_selection_scale = (1.0 + _selection_size) * exchange_scale
        exchange_selection_scaling = np.diag([exchange_selection_scale, exchange_selection_scale, 1.0, 1.0]).astype(
            np.float32
        )

        hover_set = set(self._hover_ids)
        selection_set = set(self._current_selection)

        # 1. Sites
        if self._display_options.show_sites:
            for site in self.render_model.sites:
                if not (site.is_magnetic or self._display_options.show_nonmagnetic_atoms):
                    continue

                site_color = (
                    self._display_options.default_site_color if site.color is None else site.color
                )

                if site.render_id in hover_set:
                    mode = (
                        SelectionMode.SELECTED_HOVER
                        if site.render_id in selection_set
                        else SelectionMode.HOVER
                    )
                elif site.render_id in selection_set:
                    mode = SelectionMode.SELECTED
                else:
                    mode = SelectionMode.NOT_SELECTED

                show_arrow = site.is_magnetic and self._display_options.show_spins_not_atoms
                render_object = self.arrow if show_arrow else self.small_sphere
                selection_render_object = self.selection_arrow if show_arrow else self.selection_small_sphere

                for model_matrix in site.model_matrices(self._display_options.prettify):
                    site_model_matrix = model_matrix @ spin_scale_matrix

                    if not show_arrow:
                        s = site.radius
                        atomscale = np.diag([s, s, s, 1.0]).astype(np.float32)
                        site_model_matrix = site_model_matrix @ atomscale

                    if mode != SelectionMode.NOT_SELECTED:
                        self.selection_shader.use(use_custom_proj, custom_proj)
                        self.selection_shader.set_object(site_model_matrix, mode)
                        selection_render_object.render_back(self.gl)

                    self.object_shader.use(use_custom_proj, custom_proj)
                    self.object_shader.set_object(site_model_matrix, site_color)
                    render_object.render_triangles(self.gl)

        # 2. Exchanges
        if self._display_options.show_exchanges:
            for exchange in self.render_model.exchanges:
                exchange_color = (
                    self._display_options.default_exchange_color
                    if exchange.color is None
                    else exchange.color
                )

                if exchange.render_id in hover_set:
                    mode = (
                        SelectionMode.SELECTED_HOVER
                        if exchange.render_id in selection_set
                        else SelectionMode.HOVER
                    )
                elif exchange.render_id in selection_set:
                    mode = SelectionMode.SELECTED
                else:
                    mode = SelectionMode.NOT_SELECTED

                for model_matrix in exchange.model_matrices(self._display_options.prettify):
                    if mode != SelectionMode.NOT_SELECTED:
                        exchange_sel_matrix = model_matrix @ exchange_selection_scaling
                        self.selection_shader.use(use_custom_proj, custom_proj)
                        self.selection_shader.set_object(exchange_sel_matrix, mode)
                        self.selection_tube.render_back(self.gl)

                    exchange_model_matrix = model_matrix @ exchange_scaling
                    self.object_shader.use(use_custom_proj, custom_proj)
                    self.object_shader.set_object(exchange_model_matrix, exchange_color)
                    self.tube.render_triangles(self.gl)

        # 3. Anisotropies (alpha-blended back-to-front)
        if self._display_options.show_anisotropies and self.render_model.anisotropies:
            self.gl.enable(blend=True)
            self.gl.blend_func("SRC_ALPHA", "ONE_MINUS_SRC_ALPHA")

            camera_vector = np.array(self.gl.camera_pos, dtype=float) - self.view_origin
            sorted_anisotropies = sorted(
                self.render_model.anisotropies,
                key=lambda aniso: float(np.dot(aniso.xyz, camera_vector)),
            )

            self.anisotropy_shader.use(
                self._display_options.anisotropy_color,
                use_custom_proj,
                custom_proj,
            )

            for anisotropy in sorted_anisotropies:
                for model_matrix in anisotropy.model_matrices(self._display_options.prettify):
                    anisotropy_model_matrix = model_matrix @ anisotropy_scale_matrix
                    self.anisotropy_shader.set_model_matrix(anisotropy_model_matrix)
                    self.anisotropy_sphere.render_triangles(self.gl)

            self.gl.disable(blend=True)

        # 4. Supercell Wireframe
        if self._display_options.show_supercell:
            supercell_matrix = np.eye(4, dtype=np.float32)
            supercell_matrix[:3, :3] = self.render_model.expanded.structure.unit_cell._xyz.T

            self.cell_shader.use(use_custom_proj, custom_proj)
            self.cell_shader.set_model_matrix(supercell_matrix)
            self.cube.render_wireframe(self.gl)

        # 5. Unit Cell Wireframes
        if self._display_options.show_unit_cell:
            unit_cell_transform = np.eye(4, dtype=np.float32)
            unit_cell_transform[:3, :3] = self.render_model.hamiltonian.structure.unit_cell._xyz.T

            self.cell_shader.use(use_custom_proj, custom_proj)
            for offset in self.render_model.cells():
                translation_matrix = np.eye(4, dtype=np.float32)
                translation_matrix[:3, 3] = offset.vector
                self.cell_shader.set_model_matrix(unit_cell_transform @ translation_matrix)
                self.cube.render_wireframe(self.gl)

        # 6. Bottom-left Corner Viewport: Cartesian Axes
        if self._display_options.show_cartesian_axes:
            self.gl.disable(depth_test=True)
            self.gl.viewport(
                self.axes_padding,
                self.axes_padding,
                self.axes_size,
                self.axes_size,
            )
            self.axes_shader.use()

            self.axes_shader.set_axis(_x_axis_mat, (1.0, 0.0, 0.0))
            self.arrow.render_triangles(self.gl)

            self.axes_shader.set_axis(_y_axis_mat, (0.0, 1.0, 0.0))
            self.arrow.render_triangles(self.gl)

            self.axes_shader.set_axis(_z_axis_mat, (0.0, 0.0, 1.0))
            self.arrow.render_triangles(self.gl)

            self.gl.viewport(0, 0, w, h)
            self.gl.enable(depth_test=True)

        # 7. Bottom-right Corner Viewport: Lattice Axes
        if self._display_options.show_lattice_axes:
            self.gl.disable(depth_test=True)
            self.gl.viewport(
                max(0, w - self.axes_size - self.axes_padding),
                self.axes_padding,
                self.axes_size,
                self.axes_size,
            )
            self.axes_shader.use()

            transforms = (
                self.render_model.unit_cell_orthobasis_transforms
                if self._display_options.orthogonal_lattice_axes
                else self.render_model.unit_cell_axes_transforms
            )

            self.axes_shader.set_axis(transforms[0], (0.8, 0.4, 0.0))
            self.arrow.render_triangles(self.gl)

            self.axes_shader.set_axis(transforms[1], (0.6, 0.8, 0.0))
            self.arrow.render_triangles(self.gl)

            self.axes_shader.set_axis(transforms[2], (0.7, 0.3, 0.9))
            self.arrow.render_triangles(self.gl)

            self.gl.viewport(0, 0, w, h)
            self.gl.enable(depth_test=True)

        self.gl.execute_commands(execute_once=False, clear_previous=True)

    def snapshot(self) -> np.ndarray:
        """Return the current frame as an RGBA uint8 NumPy array.

        If the frontend has synced ``image_data`` via ``sync_image_data=True``,
        returns ``self.gl.get_image_data()``. Otherwise returns a blank frame
        initialized to ``background_color`` when running headless without a browser.

        Returns
        -------
        (height, width, 4) numpy.ndarray
            RGBA ``uint8`` image array.
        """
        if self.gl.image_data is not None:
            return self.gl.get_image_data()

        h, w = int(self.gl.height), int(self.gl.width)
        frame = np.zeros((h, w, 4), dtype=np.uint8)
        r, g, b = self._display_options.background_color
        frame[:, :, 0] = int(np.clip(r * 255, 0, 255))
        frame[:, :, 1] = int(np.clip(g * 255, 0, 255))
        frame[:, :, 2] = int(np.clip(b * 255, 0, 255))
        frame[:, :, 3] = 255
        return frame
