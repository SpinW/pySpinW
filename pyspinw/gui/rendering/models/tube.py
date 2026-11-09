"""Tube Model."""

import numpy as np

from pyspinw.gui.rendering.model import Model


class Tube(Model):
    """3D Model of hollow tube.

    Parameters
    ----------
    n_points : int, optional
        Number of angular sample points around the circumference. Default is ``51``.
    padding : float, optional
        Radial padding added to the unit tube radius (e.g., for selection outlines).
        Default is ``0.0``.
    """

    def __init__(self, n_points=51, padding=0.0):
        super().__init__()

        vertices_and_normals = self.build_geometry(n_points=n_points, padding=padding)

        self.add_vertex_normal_data(vertices_and_normals)

    @staticmethod
    def build_geometry(n_points: int = 51, padding: float = 0.0) -> np.ndarray:
        """Build the vertex and normal array for a hollow tube without OpenGL calls.

        Parameters
        ----------
        n_points : int, optional
            Number of angular sample points around the circumference. Default is ``51``.
        padding : float, optional
            Radial padding added to the unit tube radius (e.g., for selection outlines).
            Default is ``0.0``.

        Returns
        -------
        (n, 6) numpy.ndarray
            2D ``float32`` array of interleaved 3D vertex positions and 3D normal vectors.
        """
        angles = np.linspace(0, 2*np.pi, n_points)

        vertices_and_normals = []

        r = 1.0 + padding

        for i in range(n_points-1):
            angle_0 = angles[i]
            angle_1 = angles[i+1]

            x_0 = np.cos(angle_0)
            x_1 = np.cos(angle_1)

            y_0 = np.sin(angle_0)
            y_1 = np.sin(angle_1)

            vertices_and_normals += [
                [x_0*r, y_0*r, 0, x_0, y_0, 0],
                [x_1*r, y_1*r, 0, x_1, y_1, 0],
                [x_0*r, y_0*r, 1, x_0, y_0, 0],
                [x_1*r, y_1*r, 0, x_1, y_1, 0],
                [x_1*r, y_1*r, 1, x_1, y_1, 0],
                [x_0*r, y_0*r, 1, x_0, y_0, 0],
            ]

        return np.array(vertices_and_normals, dtype=np.float32)
