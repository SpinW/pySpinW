import numpy as np

from pyspinw.gui.rendering.model import Model


class Arrow(Model):
    def __init__(self, n_points=51):
        super().__init__()

        angles = np.linspace(0, 2*np.pi, n_points)

        vertices_and_normals = []

        for i in range(n_points-1):
            angle_0 = angles[i]
            angle_1 = angles[i+1]

            x_0 = np.cos(angle_0)
            x_1 = np.cos(angle_1)

            y_0 = np.sin(angle_0)
            y_1 = np.sin(angle_1)

            vertices_and_normals += [
                [x_0, y_0, 0, x_0, y_0, 0],
                [x_1, y_1, 0, x_1, y_1, 0],
                [x_0, y_0, 1, x_0, y_0, 0],
                [x_1, y_1, 0, x_1, y_1, 0],
                [x_1, y_1, 1, x_1, y_1, 0],
                [x_0, y_0, 1, x_0, y_0, 0],
            ]

        vertices_and_normals = np.array(vertices_and_normals, dtype=np.float32)

        self.add_vertex_normal_data(vertices_and_normals)