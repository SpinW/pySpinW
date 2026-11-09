"""Sphere Model."""

import numpy as np
from scipy.spatial import ConvexHull

from pyspinw.calculations.geodesic import Geodesic
from pyspinw.gui.rendering.model import Model


class Sphere(Model):
    """3D Model of a sphere.

    Parameters
    ----------
    divisions : int
        Number of geodesic subdivisions used to generate the sphere mesh.
    radius : float, optional
        Radius of the sphere. Default is ``1.0``.
    padding : float, optional
        Radial padding added to the sphere radius (e.g., for selection outlines).
        Default is ``0.0``.

    Attributes
    ----------
    vertices_and_normals : numpy.ndarray
        Flattened ``float32`` array of interleaved vertex positions and normals.
    """

    def __init__(self, divisions: int, radius: float=1.0, padding=0.0):

        super().__init__()

        self.vertices_and_normals = self.build_geometry(
            divisions=divisions,
            radius=radius,
            padding=padding,
        ).reshape(-1)

        self.add_vertex_normal_data(self.vertices_and_normals)

    @staticmethod
    def build_geometry(divisions: int, radius: float = 1.0, padding: float = 0.0) -> np.ndarray:
        """Build the vertex and normal array for a geodesic sphere without OpenGL calls.

        Parameters
        ----------
        divisions : int
            Number of geodesic subdivisions used to generate the sphere mesh.
        radius : float, optional
            Radius of the sphere. Default is ``1.0``.
        padding : float, optional
            Radial padding added to the sphere radius (e.g., for selection outlines).
            Default is ``0.0``.

        Returns
        -------
        (n, 6) numpy.ndarray
            2D ``float32`` array of interleaved 3D vertex positions and 3D normal vectors.
        """
        # Build list of vertices
        points, _ = Geodesic.by_divisions(divisions)

        # Need to arrange these into faces
        hull = ConvexHull(points)

        faces = hull.simplices

        verts = []
        for raw_face in faces:

            # Dereference faces
            face = hull.vertices[raw_face]

            # Find orientation
            a = points[face[1], :] - points[face[0], :]
            b = points[face[2], :] - points[face[0], :]
            n = np.cross(a, b)

            c = (points[face[0], :] + points[face[1], :] + points[face[2], :]) / 3

            if np.dot(n, c) < 0:
                face = face[::-1]

            # Add verts
            for i in face:
                verts.append(np.concatenate(((radius + padding) * points[i, :], points[i, :])))

        return np.array(verts, dtype=np.float32)


if __name__ == "__main__":
    s = Sphere(2)
