import numpy as np
from scipy.spatial import ConvexHull

from pyspinw.calculations.geodesic import Geodesic


class Sphere:
    def __init__(self, divisions: int):
        self._divisions = divisions

        self.vertices_and_normals = np.array([[]], dtype=float)

        self._recalculate_vertices()



    def _recalculate_vertices(self):
        # Build list of vertices

        points, _ = Geodesic.by_divisions(self._divisions)

        # Need to arrange these into faces

        hull = ConvexHull(points)

        faces = hull.simplices

        verts = []
        for face in faces:

            # Dereference faces
            face = hull.vertices[face]


            # Find orientation
            a = points[face[1], :] - points[face[0], :]
            b = points[face[2], :] - points[face[0], :]
            n = np.cross(a, b)

            c = (points[face[0], :] + points[face[1], :] + points[face[2], :]) / 3

            if np.dot(n, c) > 0:
                face = face[::-1]

            # Add verts
            for i in face:
                verts.append(points[i, :])
                verts.append(points[i, :]/np.sqrt(np.sum(points[i, :]**2))) # Normals, normalised


        self.vertices_and_normals = np.array(verts, dtype=np.float32).reshape(-1) # Return as float32 for rendering

    @property
    def divisions(self):
        return self._divisions

    @divisions.setter
    def divisions(self, divisions):
        self._divisions = divisions
        self._recalculate_vertices()


if __name__ == "__main__":
    s = Sphere(2)