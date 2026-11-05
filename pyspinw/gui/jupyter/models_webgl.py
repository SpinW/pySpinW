"""Mesh generation and WebGL 2.0 buffer/VAO management for Jupyter notebooks."""

from __future__ import annotations

import numpy as np
from scipy.spatial import ConvexHull

from pyspinw.calculations.geodesic import Geodesic


# ---------------------------------------------------------------------------
# Pure NumPy / Scipy Geometry Generators (Free of PyOpenGL / Qt)
# ---------------------------------------------------------------------------

def generate_sphere_geometry(divisions: int, radius: float = 1.0, padding: float = 0.0) -> np.ndarray:
    """Generate vertex and normal data for a geodesic sphere.

    Parameters
    ----------
    divisions : int
        Number of geodesic subdivisions.
    radius : float, optional
        Radius of the sphere. Default is 1.0.
    padding : float, optional
        Radial padding added to the sphere radius. Default is 0.0.

    Returns
    -------
    numpy.ndarray
        Array of shape (N, 6) with float32 vertex coordinates and normals.
    """
    points, _ = Geodesic.by_divisions(divisions)
    hull = ConvexHull(points)

    faces = hull.simplices
    verts = []
    for raw_face in faces:
        face = hull.vertices[raw_face]

        # Orientation check
        a = points[face[1], :] - points[face[0], :]
        b = points[face[2], :] - points[face[0], :]
        n = np.cross(a, b)
        c = (points[face[0], :] + points[face[1], :] + points[face[2], :]) / 3.0

        if np.dot(n, c) < 0:
            face = face[::-1]

        for i in face:
            verts.append((radius + padding) * points[i, :])
            verts.append(points[i, :])  # Normal

    arr = np.array(verts, dtype=np.float32).reshape(-1, 6)
    return arr


def generate_tube_geometry(n_points: int = 51, padding: float = 0.0) -> np.ndarray:
    """Generate vertex and normal data for a cylindrical tube.

    Parameters
    ----------
    n_points : int, optional
        Number of circumference sample points. Default is 51.
    padding : float, optional
        Radial padding added to the tube radius. Default is 0.0.

    Returns
    -------
    numpy.ndarray
        Array of shape (N, 6) with float32 vertex coordinates and normals.
    """
    angles = np.linspace(0, 2 * np.pi, n_points)
    r = 1.0 + padding
    verts = []

    for i in range(n_points - 1):
        a0 = angles[i]
        a1 = angles[i + 1]

        x0 = np.cos(a0)
        x1 = np.cos(a1)
        y0 = np.sin(a0)
        y1 = np.sin(a1)

        verts += [
            [x0 * r, y0 * r, 0.0, x0, y0, 0.0],
            [x1 * r, y1 * r, 0.0, x1, y1, 0.0],
            [x0 * r, y0 * r, 1.0, x0, y0, 0.0],
            [x1 * r, y1 * r, 0.0, x1, y1, 0.0],
            [x1 * r, y1 * r, 1.0, x1, y1, 0.0],
            [x0 * r, y0 * r, 1.0, x0, y0, 0.0],
        ]

    return np.array(verts, dtype=np.float32)


def generate_arrow_geometry(
    n_points: int = 51,
    r_shaft: float = 0.05,
    r_head: float = 0.1,
    head_fraction: float = 0.3,
    padding: float = 0.0,
    modify_ring: bool = False,
) -> np.ndarray:
    """Generate vertex and normal data for a 3D arrow centered at the origin.

    Parameters
    ----------
    n_points : int, optional
        Number of angular samples. Default is 51.
    r_shaft : float, optional
        Radius of shaft cylinder. Default is 0.05.
    r_head : float, optional
        Radius of cone base. Default is 0.1.
    head_fraction : float, optional
        Fraction of length taken by cone head. Default is 0.3.
    padding : float, optional
        Radial/axial padding. Default is 0.0.
    modify_ring : bool, optional
        Whether to flip ring normals radially inward for outline culling. Default is False.

    Returns
    -------
    numpy.ndarray
        Array of shape (N, 6) with float32 vertex coordinates and normals.
    """
    angles = np.linspace(0, 2 * np.pi, n_points)
    mid_angles = 0.5 * (angles[1:] + angles[:-1])

    x = np.cos(angles)
    y = np.sin(angles)

    x_head = (r_head + padding) * x
    y_head = (r_head + padding) * y

    x_shaft = (r_shaft + padding) * x
    y_shaft = (r_shaft + padding) * y

    z_transition = 1.0 - (head_fraction + 0.01 * padding)
    verts = []

    # 1. Cone Head
    rz = np.array([r_head, head_fraction])
    head_normal_z, head_normal_r = rz / np.sqrt(np.sum(rz**2))
    head_normal_z *= -1.0

    x_mid = head_normal_r * np.cos(mid_angles)
    y_mid = head_normal_r * np.sin(mid_angles)
    z_with_pad = 1.0 + padding * head_fraction / r_head

    for i in range(n_points - 1):
        x0 = x_head[i]
        x1 = x_head[i + 1]
        y0 = y_head[i]
        y1 = y_head[i + 1]

        nx0 = head_normal_r * x[i]
        nx1 = head_normal_r * x[i + 1]
        ny0 = head_normal_r * y[i]
        ny1 = head_normal_r * y[i + 1]

        verts += [
            [x0, y0, z_transition, nx0, ny0, head_normal_z],
            [x1, y1, z_transition, nx1, ny1, head_normal_z],
            [0.0, 0.0, z_with_pad, x_mid[i], y_mid[i], head_normal_z],
        ]

    # 2. Ring base between shaft and head
    for i in range(n_points - 1):
        x_in0, x_in1 = x_shaft[i], x_shaft[i + 1]
        y_in0, y_in1 = y_shaft[i], y_shaft[i + 1]
        x_out0, x_out1 = x_head[i], x_head[i + 1]
        y_out0, y_out1 = y_head[i], y_head[i + 1]

        if modify_ring:
            nx0, nx1 = -x[i], -x[i + 1]
            ny0, ny1 = -y[i], -y[i + 1]
            nz = 0.0
        else:
            nx0, nx1 = 0.0, 0.0
            ny0, ny1 = 0.0, 0.0
            nz = -1.0

        verts += [
            [x_in0, y_in0, z_transition, nx0, ny0, nz],
            [x_in1, y_in1, z_transition, nx1, ny1, nz],
            [x_out0, y_out0, z_transition, nx0, ny0, nz],
            [x_in1, y_in1, z_transition, nx1, ny1, nz],
            [x_out1, y_out1, z_transition, nx1, ny1, nz],
            [x_out0, y_out0, z_transition, nx0, ny0, nz],
        ]

    # 3. Shaft cylinder
    for i in range(n_points - 1):
        x0, x1 = x_shaft[i], x_shaft[i + 1]
        y0, y1 = y_shaft[i], y_shaft[i + 1]
        nx0, nx1 = x[i], x[i + 1]
        ny0, ny1 = y[i], y[i + 1]

        verts += [
            [x0, y0, 0.0, nx0, ny0, 0.0],
            [x1, y1, 0.0, nx1, ny1, 0.0],
            [x0, y0, z_transition, nx0, ny0, 0.0],
            [x1, y1, 0.0, nx1, ny1, 0.0],
            [x1, y1, z_transition, nx1, ny1, 0.0],
            [x0, y0, z_transition, nx0, ny0, 0.0],
        ]

    # 4. Base cap
    for i in range(n_points - 1):
        x0, x1 = x_shaft[i], x_shaft[i + 1]
        y0, y1 = y_shaft[i], y_shaft[i + 1]

        verts += [
            [x1, y1, 0.0, 0.0, 0.0, -1.0],
            [x0, y0, 0.0, 0.0, 0.0, -1.0],
            [-padding, -padding, 0.0, 0.0, 0.0, -1.0],
        ]

    arr = np.array(verts, dtype=np.float32)

    # Center arrow at the origin along z: [-0.5, 0.5]
    arr[:, 2] -= 0.5
    return arr


def generate_cube_wireframe_geometry() -> np.ndarray:
    """Generate 24 line vertices for a unit cube wireframe [0, 1]^3.

    Returns
    -------
    numpy.ndarray
        Array of shape (24, 3) with float32 3D line coordinates.
    """
    edges = np.array(
        [
            [0, 0, 0], [0, 0, 1],  # Z direction
            [1, 0, 0], [1, 0, 1],
            [0, 1, 0], [0, 1, 1],
            [1, 1, 0], [1, 1, 1],
            [0, 0, 0], [0, 1, 0],  # Y direction
            [1, 0, 0], [1, 1, 0],
            [0, 0, 1], [0, 1, 1],
            [1, 0, 1], [1, 1, 1],
            [0, 0, 0], [1, 0, 0],  # X direction
            [0, 1, 0], [1, 1, 0],
            [0, 0, 1], [1, 0, 1],
            [0, 1, 1], [1, 1, 1],
        ],
        dtype=np.float32,
    )
    return edges


# ---------------------------------------------------------------------------
# WebGL Mesh / Buffer Class
# ---------------------------------------------------------------------------

class WebGLMesh:
    """Wrapper that binds vertex data to a WebGL2 buffer and VAO on ipywebgl.GLViewer.

    Parameters
    ----------
    viewer : ipywebgl.GLViewer
        The WebGL viewer widget.
    vertex_data : numpy.ndarray
        Vertex data array (N, 6) for position+normal or (N, 3) for position only.
    has_normals : bool, optional
        Whether normal vectors are present (stride = 6 * 4). Default is True.
    draw_mode : str, optional
        WebGL draw primitive ('TRIANGLES' or 'LINES'). Default is 'TRIANGLES'.
    """

    def __init__(
        self,
        viewer,
        vertex_data: np.ndarray,
        has_normals: bool = True,
        draw_mode: str = "TRIANGLES",
    ):
        self.viewer = viewer
        self.draw_mode = draw_mode
        self.has_normals = has_normals

        data = np.ascontiguousarray(vertex_data, dtype=np.float32)
        if has_normals:
            self.vertex_count = len(data.reshape(-1, 6))
            stride = 6 * 4
        else:
            self.vertex_count = len(data.reshape(-1, 3))
            stride = 3 * 4

        self.vbo = viewer.create_buffer()
        viewer.bind_buffer("ARRAY_BUFFER", self.vbo)
        viewer.buffer_data("ARRAY_BUFFER", data, "STATIC_DRAW")

        self.vao = viewer.create_vertex_array()
        viewer.bind_vertex_array(self.vao)

        # Attribute 0: position
        viewer.enable_vertex_attrib_array(0)
        viewer.vertex_attrib_pointer(0, 3, "FLOAT", False, stride, 0)

        # Attribute 1: normal (if present)
        if has_normals:
            viewer.enable_vertex_attrib_array(1)
            viewer.vertex_attrib_pointer(1, 3, "FLOAT", False, stride, 3 * 4)

        viewer.bind_vertex_array(None)
        viewer.bind_buffer("ARRAY_BUFFER", None)

    def render_triangles(self):
        """Render using normal back-face culling."""
        self.viewer.enable("CULL_FACE")
        self.viewer.cull_face("BACK")
        self.viewer.bind_vertex_array(self.vao)
        self.viewer.draw_arrays(self.draw_mode, 0, self.vertex_count)
        self.viewer.bind_vertex_array(None)

    def render_back(self):
        """Render back faces only (used for selection highlighting outlines)."""
        self.viewer.enable("CULL_FACE")
        self.viewer.cull_face("FRONT")
        self.viewer.bind_vertex_array(self.vao)
        self.viewer.draw_arrays(self.draw_mode, 0, self.vertex_count)
        self.viewer.bind_vertex_array(None)

    def render_wireframe(self):
        """Render wireframe lines."""
        self.viewer.bind_vertex_array(self.vao)
        self.viewer.draw_arrays(self.draw_mode, 0, self.vertex_count)
        self.viewer.bind_vertex_array(None)

