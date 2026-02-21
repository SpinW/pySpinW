""" Sphere Model """

import numpy as np
from scipy.spatial import ConvexHull

from pyspinw.calculations.geodesic import Geodesic
from pyspinw.gui.rendering.model import Model


class Sphere(Model):
    """ 3D Model of a sphere """

    def __init__(self, divisions: int, radius: float=1.0):

        super().__init__()

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
                verts.append(radius * points[i, :])
                verts.append(points[i, :]) # Normals


        self.vertices_and_normals = np.array(verts, dtype=np.float32).reshape(-1) # Return as float32 for rendering

        print(self.vertices_and_normals.reshape(-1, 6))

        self.add_vertex_normal_data(self.vertices_and_normals)


if __name__ == "__main__":
    s = Sphere(2)
