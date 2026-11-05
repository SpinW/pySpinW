"""Jupyter Notebook widget for rendering magnetic crystal structures using ipywidgets and ipywebgl."""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING
import ipywidgets as widgets
import numpy as np

from pyspinw.gui.camera import Camera
from pyspinw.gui.displayoptions import DisplayOptions
from pyspinw.gui.jupyter import check_jupyter_dependencies
from pyspinw.gui.jupyter.interaction import CanvasInteractionHandler
from pyspinw.gui.jupyter.models_webgl import (
    WebGLMesh,
    generate_arrow_geometry,
    generate_cube_wireframe_geometry,
    generate_sphere_geometry,
    generate_tube_geometry,
)
from pyspinw.gui.jupyter.shaders_webgl import (
    AnisotropyShaderWebGL,
    AxesShaderWebGL,
    CellShaderWebGL,
    ObjectShaderWebGL,
    SelectionShaderWebGL,
)
from pyspinw.gui.rendering.selectionmode import SelectionMode
from pyspinw.util import rotation_from_z

if TYPE_CHECKING:
    from pyspinw.gui.rendermodel import RenderModel

logger = logging.getLogger(__name__)

# Transforms for xyz arrows in axes viewport
_Z_AXIS_MAT = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0.5], [0, 0, 0, 1]], dtype=np.float32)
_Y_AXIS_MAT = np.array([[1, 0, 0, 0], [0, 0, 1, 0.5], [0, -1, 0, 0], [0, 0, 0, 1]], dtype=np.float32)
_X_AXIS_MAT = np.array([[0, 0, 1, 0.5], [0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1]], dtype=np.float32)

_SELECTION_SIZE = 0.2


class CrystalViewerWidget(widgets.VBox):
    """Jupyter widget to show and interact with magnetic crystal structures.

    Parameters
    ----------
    render_model : RenderModel
        Model containing the crystal structure and magnetic interactions to render.
    initial_rotation : (3, 3) numpy.ndarray or None, optional
        Initial 3x3 view rotation matrix. Default is None.
    initial_distance : float or None, optional
        Initial camera distance from the view origin. Default is None (10.0).
    display_options : DisplayOptions or None, optional
        Initial display options. Default is None.
    width : int, optional
        Viewport width in pixels. Default is 700.
    height : int, optional
        Viewport height in pixels. Default is 500.
    """

    def __init__(
        self,
        render_model: RenderModel,
        initial_rotation: np.ndarray | None = None,
        initial_distance: float | None = None,
        display_options: DisplayOptions | None = None,
        width: int = 700,
        height: int = 500,
    ):
        check_jupyter_dependencies()
        super().__init__()

        import ipywebgl

        self.width = width
        self.height = height
        self.render_model = render_model
        self.display_options = display_options if display_options is not None else DisplayOptions()
        self.camera = Camera(horizontal_pixels=width, vertical_pixels=height)

        # View details
        self.view_origin = np.array(render_model.expanded.structure.unit_cell.centre, dtype=float)
        self.view_rotation = np.eye(3, dtype=float) if initial_rotation is None else np.array(initial_rotation, dtype=float)
        self.view_radius = 10.0 if initial_distance is None else float(initial_distance)

        # Mouse transformation deltas
        self.mouse_rotation = np.eye(3, dtype=float)
        self.mouse_origin = np.zeros(3, dtype=float)

        # Selection state
        self.hover_ids: list[int] = []
        self.current_selection: list[int] = []

        # Create WebGL viewer canvas
        self.gl_viewer = ipywebgl.GLViewer(width=width, height=height)

        # Attach mouse interaction handler
        self.interaction = CanvasInteractionHandler(self)

        # Setup WebGL shaders and meshes
        self._init_gl()

        # Set widget children
        self.children = [self.gl_viewer]

        # Render the initial frame
        self.render()

    def _init_gl(self) -> None:
        """Compile shaders and create geometry meshes on the WebGL viewer."""
        # Shaders
        self.object_shader = ObjectShaderWebGL(self.gl_viewer)
        self.selection_shader = SelectionShaderWebGL(self.gl_viewer)
        self.cell_shader = CellShaderWebGL(self.gl_viewer)
        self.anisotropy_shader = AnisotropyShaderWebGL(self.gl_viewer)
        self.axes_shader = AxesShaderWebGL(self.gl_viewer)

        # Meshes
        self.sphere = WebGLMesh(self.gl_viewer, generate_sphere_geometry(3))
        self.selection_sphere = WebGLMesh(self.gl_viewer, generate_sphere_geometry(3, padding=0.1))

        self.small_sphere = WebGLMesh(self.gl_viewer, generate_sphere_geometry(3, radius=0.1))
        self.selection_small_sphere = WebGLMesh(
            self.gl_viewer, generate_sphere_geometry(3, radius=0.1, padding=0.01)
        )
        self.anisotropy_sphere = WebGLMesh(self.gl_viewer, generate_sphere_geometry(5, radius=0.1))

        self.tube = WebGLMesh(self.gl_viewer, generate_tube_geometry())
        self.selection_tube = WebGLMesh(self.gl_viewer, generate_tube_geometry(padding=0.2))

        self.arrow = WebGLMesh(self.gl_viewer, generate_arrow_geometry())
        self.selection_arrow = WebGLMesh(
            self.gl_viewer, generate_arrow_geometry(padding=0.03, modify_ring=True)
        )

        self.cube = WebGLMesh(
            self.gl_viewer,
            generate_cube_wireframe_geometry(),
            has_normals=False,
            draw_mode="LINES",
        )

    def request_render(self) -> None:
        """Request a redraw of the scene."""
        self.render()

    def reset_view(self) -> None:
        """Reset view to default orientation and origin."""
        self.view_origin = np.array(self.render_model.expanded.structure.unit_cell.centre, dtype=float)
        self.view_rotation = np.eye(3, dtype=float)
        self.view_radius = 10.0
        self.mouse_rotation = np.eye(3, dtype=float)
        self.mouse_origin = np.zeros(3, dtype=float)
        self.render()

    def align_view(self, axis: str) -> None:
        """Align view direction along a crystallographic or Cartesian axis.

        Parameters
        ----------
        axis : str
            One of 'a', 'b', 'c', 'a*', 'b*', 'c*', 'x', 'y', 'z', '-x', '-y', '-z'.
        """
        unit_cell = self.render_model.hamiltonian.structure.unit_cell
        axis_lower = axis.strip().lower()

        if axis_lower in ("x", "+x"):
            v = np.array([1.0, 0.0, 0.0])
        elif axis_lower == "-x":
            v = np.array([-1.0, 0.0, 0.0])
        elif axis_lower in ("y", "+y"):
            v = np.array([0.0, 1.0, 0.0])
        elif axis_lower == "-y":
            v = np.array([0.0, -1.0, 0.0])
        elif axis_lower in ("z", "+z"):
            v = np.array([0.0, 0.0, 1.0])
        elif axis_lower == "-z":
            v = np.array([0.0, 0.0, -1.0])
        elif axis_lower in ("a", "+a"):
            v = np.array(unit_cell.lattice_units_to_cartesian([1.0, 0.0, 0.0]), dtype=float)
        elif axis_lower == "-a":
            v = -np.array(unit_cell.lattice_units_to_cartesian([1.0, 0.0, 0.0]), dtype=float)
        elif axis_lower in ("b", "+b"):
            v = np.array(unit_cell.lattice_units_to_cartesian([0.0, 1.0, 0.0]), dtype=float)
        elif axis_lower == "-b":
            v = -np.array(unit_cell.lattice_units_to_cartesian([0.0, 1.0, 0.0]), dtype=float)
        elif axis_lower in ("c", "+c"):
            v = np.array(unit_cell.lattice_units_to_cartesian([0.0, 0.0, 1.0]), dtype=float)
        elif axis_lower == "-c":
            v = -np.array(unit_cell.lattice_units_to_cartesian([0.0, 0.0, 1.0]), dtype=float)
        elif axis_lower in ("a*", "+a*"):
            b_cart = np.array(unit_cell.lattice_units_to_cartesian([0.0, 1.0, 0.0]), dtype=float)
            c_cart = np.array(unit_cell.lattice_units_to_cartesian([0.0, 0.0, 1.0]), dtype=float)
            v = np.cross(b_cart, c_cart)
        elif axis_lower in ("b*", "+b*"):
            a_cart = np.array(unit_cell.lattice_units_to_cartesian([1.0, 0.0, 0.0]), dtype=float)
            c_cart = np.array(unit_cell.lattice_units_to_cartesian([0.0, 0.0, 1.0]), dtype=float)
            v = np.cross(c_cart, a_cart)
        elif axis_lower in ("c*", "+c*"):
            a_cart = np.array(unit_cell.lattice_units_to_cartesian([1.0, 0.0, 0.0]), dtype=float)
            b_cart = np.array(unit_cell.lattice_units_to_cartesian([0.0, 1.0, 0.0]), dtype=float)
            v = np.cross(a_cart, b_cart)
        else:
            raise ValueError(f"Unknown axis: {axis}. Supported: 'a', 'b', 'c', 'a*', 'b*', 'c*', 'x', 'y', 'z'.")

        norm = np.linalg.norm(v)
        if norm > 1e-8:
            v = v / norm

        self.view_rotation = rotation_from_z(v)
        self.mouse_rotation = np.eye(3, dtype=float)
        self.mouse_origin = np.zeros(3, dtype=float)
        self.render()

    def snapshot(self) -> np.ndarray | None:
        """Capture the current WebGL canvas buffer if available."""
        if hasattr(self.gl_viewer, "image_data"):
            return self.gl_viewer.image_data
        return None

    def on_canvas_click(self, shift: bool = False) -> None:
        """Handle canvas click for element selection."""
        # Can be connected to element picking or selection callbacks
        pass

    def render(self) -> None:
        """Render the complete crystal structure and magnetic interactions."""
        # 1. Update camera coordinates
        compound_rotation = self.view_rotation @ self.mouse_rotation
        origin_world = self.view_origin + self.mouse_origin
        camera_world = origin_world + compound_rotation @ np.array([0, 0, self.view_radius], dtype=float)
        up_world = compound_rotation @ np.array([0, 1, 0], dtype=float)

        self.camera.position = tuple(camera_world)
        self.camera.up = tuple(up_world)
        self.camera.look_at = tuple(origin_world)
        self.camera.orthographic = not self.display_options.perspective
        self.camera.fixed_size = None

        # 2. Viewport & Buffer Clear
        pw, ph = self.width, self.height
        self.gl_viewer.viewport(0, 0, pw, ph)
        self.gl_viewer.enable("DEPTH_TEST")
        self.gl_viewer.enable("CULL_FACE")
        self.gl_viewer.disable("BLEND")

        bg_r, bg_g, bg_b = self.display_options.background_color
        self.gl_viewer.clear_color(bg_r, bg_g, bg_b, 1.0)
        self.gl_viewer.clear()

        # Scaling matrices
        spin_scale = 2.0 * self.display_options.atom_spin_scaling
        spin_scale_matrix = np.diag([spin_scale, spin_scale, spin_scale, 1.0])

        anisotropy_scale = 3.0 * self.display_options.anisotropy_scaling
        anisotropy_scale_matrix = np.diag([anisotropy_scale, anisotropy_scale, anisotropy_scale, 1.0])

        exchange_scale = 0.1 * self.display_options.exchange_scaling
        exchange_scaling = np.diag([exchange_scale, exchange_scale, 1.0, 1.0])
        exchange_selection_scale = ((1.0 + _SELECTION_SIZE) / 1.0) * exchange_scale
        exchange_selection_scaling = np.diag([exchange_selection_scale, exchange_selection_scale, 1.0, 1.0])

        # 3. Render Sites (Atoms & Magnetic Spins)
        if self.display_options.show_sites:
            self.object_shader.use()
            self.object_shader.set_camera(self.camera)

            for site in self.render_model.sites:
                site_color = (
                    self.display_options.default_site_color if site.color is None else site.color
                )
                self.object_shader.set_color(site_color)

                if site.is_magnetic or self.display_options.show_nonmagnetic_atoms:
                    if site.render_id in self.hover_ids:
                        mode = (
                            SelectionMode.SELECTED_HOVER
                            if site.render_id in self.current_selection
                            else SelectionMode.HOVER
                        )
                    else:
                        mode = (
                            SelectionMode.SELECTED
                            if site.render_id in self.current_selection
                            else SelectionMode.NOT_SELECTED
                        )

                    show_arrow = site.is_magnetic and self.display_options.show_spins_not_atoms
                    render_mesh = self.arrow if show_arrow else self.small_sphere
                    selection_mesh = self.selection_arrow if show_arrow else self.selection_small_sphere

                    for model_matrix in site.model_matrices(self.display_options.prettify):
                        site_model_matrix = model_matrix @ spin_scale_matrix

                        if not show_arrow:
                            s = site.radius
                            site_model_matrix = site_model_matrix @ np.diag([s, s, s, 1.0])

                        if mode != SelectionMode.NOT_SELECTED:
                            self.selection_shader.use()
                            self.selection_shader.set_camera(self.camera)
                            self.selection_shader.set_model_matrix(site_model_matrix)
                            if mode == SelectionMode.SELECTED:
                                self.selection_shader.set_color(self.display_options.selected_color)
                            elif mode == SelectionMode.HOVER:
                                self.selection_shader.set_color(self.display_options.hover_color)
                            else:
                                self.selection_shader.set_color(self.display_options.selected_hover_color)
                            selection_mesh.render_back()

                        self.object_shader.use()
                        self.object_shader.set_model_matrix(site_model_matrix)
                        render_mesh.render_triangles()

        # 4. Render Exchanges (Bonds)
        if self.display_options.show_exchanges:
            self.object_shader.use()
            self.object_shader.set_camera(self.camera)

            for exchange in self.render_model.exchanges:
                exchange_color = (
                    self.display_options.default_exchange_color
                    if exchange.color is None
                    else exchange.color
                )
                self.object_shader.set_color(exchange_color)

                if exchange.render_id in self.hover_ids:
                    mode = (
                        SelectionMode.SELECTED_HOVER
                        if exchange.render_id in self.current_selection
                        else SelectionMode.HOVER
                    )
                else:
                    mode = (
                        SelectionMode.SELECTED
                        if exchange.render_id in self.current_selection
                        else SelectionMode.NOT_SELECTED
                    )

                for model_matrix in exchange.model_matrices(self.display_options.prettify):
                    if mode != SelectionMode.NOT_SELECTED:
                        self.selection_shader.use()
                        self.selection_shader.set_camera(self.camera)
                        self.selection_shader.set_model_matrix(model_matrix @ exchange_selection_scaling)
                        if mode == SelectionMode.SELECTED:
                            self.selection_shader.set_color(self.display_options.selected_color)
                        elif mode == SelectionMode.HOVER:
                            self.selection_shader.set_color(self.display_options.hover_color)
                        else:
                            self.selection_shader.set_color(self.display_options.selected_hover_color)
                        self.selection_tube.render_back()

                    self.object_shader.use()
                    self.object_shader.set_model_matrix(model_matrix @ exchange_scaling)
                    self.tube.render_triangles()

        # 5. Render Anisotropies
        if self.display_options.show_anisotropies and len(self.render_model.anisotropies) > 0:
            self.gl_viewer.enable("BLEND")
            self.gl_viewer.blend_func("SRC_ALPHA", "ONE_MINUS_SRC_ALPHA")

            self.anisotropy_shader.use()
            self.anisotropy_shader.set_camera(self.camera)
            self.anisotropy_shader.set_color(self.display_options.anisotropy_color)

            camera_vector = np.array(self.camera.position) - np.array(self.camera.look_at)
            sorted_anisotropies = sorted(
                self.render_model.anisotropies,
                key=lambda aniso: np.dot(aniso.xyz, camera_vector),
            )

            for anisotropy in sorted_anisotropies:
                for model_matrix in anisotropy.model_matrices(self.display_options.prettify):
                    anisotropy_model_matrix = model_matrix @ anisotropy_scale_matrix
                    self.anisotropy_shader.set_model_matrix(anisotropy_model_matrix)
                    self.anisotropy_sphere.render_triangles()

            self.gl_viewer.disable("BLEND")

        # 6. Render Supercell Wireframe
        if self.display_options.show_supercell:
            supercell_matrix = np.eye(4, dtype=np.float32)
            supercell_matrix[:3, :3] = self.render_model.expanded.structure.unit_cell._xyz.T

            self.cell_shader.use()
            self.cell_shader.set_camera(self.camera)
            self.cell_shader.set_model_matrix(supercell_matrix)
            self.cell_shader.set_color((0.7, 0.7, 0.7))
            self.cube.render_wireframe()

        # 7. Render Unit Cell Wireframes
        if self.display_options.show_unit_cell:
            unit_cell_transform = np.eye(4, dtype=np.float32)
            unit_cell_transform[:3, :3] = self.render_model.hamiltonian.structure.unit_cell._xyz.T

            self.cell_shader.use()
            self.cell_shader.set_camera(self.camera)
            self.cell_shader.set_color((0.5, 0.5, 0.5))

            for offset in self.render_model.cells():
                translation_matrix = np.eye(4, dtype=np.float32)
                translation_matrix[:3, 3] = offset.vector
                self.cell_shader.set_model_matrix(unit_cell_transform @ translation_matrix)
                self.cube.render_wireframe()

        # 8. Corner Viewports: Cartesian & Lattice Axes
        axes_size = 100
        axes_padding = 10

        if self.display_options.show_cartesian_axes:
            self.gl_viewer.disable("DEPTH_TEST")
            self.gl_viewer.viewport(axes_padding, axes_padding, axes_size, axes_size)

            self.axes_shader.use()
            self.axes_shader.set_camera(self.camera)

            # X: Red
            self.axes_shader.set_color((1.0, 0.0, 0.0))
            self.axes_shader.set_model_matrix(_X_AXIS_MAT)
            self.arrow.render_triangles()

            # Y: Green
            self.axes_shader.set_color((0.0, 1.0, 0.0))
            self.axes_shader.set_model_matrix(_Y_AXIS_MAT)
            self.arrow.render_triangles()

            # Z: Blue
            self.axes_shader.set_color((0.0, 0.0, 1.0))
            self.axes_shader.set_model_matrix(_Z_AXIS_MAT)
            self.arrow.render_triangles()

            self.gl_viewer.enable("DEPTH_TEST")
            self.gl_viewer.viewport(0, 0, pw, ph)

        if self.display_options.show_lattice_axes:
            self.gl_viewer.disable("DEPTH_TEST")
            self.gl_viewer.viewport(pw - axes_size - axes_padding, axes_padding, axes_size, axes_size)

            self.axes_shader.use()
            self.axes_shader.set_camera(self.camera)

            if self.display_options.orthogonal_lattice_axes:
                transforms = self.render_model.unit_cell_orthobasis_transforms
            else:
                transforms = self.render_model.unit_cell_axes_transforms

            # a: Orange
            self.axes_shader.set_color((0.8, 0.4, 0.0))
            self.axes_shader.set_model_matrix(transforms[0])
            self.arrow.render_triangles()

            # b: Yellow-green
            self.axes_shader.set_color((0.6, 0.8, 0.0))
            self.axes_shader.set_model_matrix(transforms[1])
            self.arrow.render_triangles()

            # c: Purple
            self.axes_shader.set_color((0.7, 0.3, 0.9))
            self.axes_shader.set_model_matrix(transforms[2])
            self.arrow.render_triangles()

            self.gl_viewer.enable("DEPTH_TEST")
            self.gl_viewer.viewport(0, 0, pw, ph)

        # 9. Flush commands to WebGL context
        self.gl_viewer.execute_commands(execute_once=True)

