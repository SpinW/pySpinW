""" Arrow Model """

import numpy as np

from pyspinw.gui.rendering.model import Model


class Arrow(Model):
    """ Model of an arrow """

    def __init__(self,
                 n_points=51,
                 r_shaft: float = 0.05,
                 r_head: float = 0.1,
                 head_fraction: float = 0.3):

        super().__init__()



        angles = np.linspace(0, 2*np.pi, n_points)

        mid_angles = 0.5*(angles[1:] + angles[:-1])

        x = np.cos(angles)
        y = np.sin(angles)

        x_head = r_head * x
        y_head = r_head * y

        x_shaft = r_shaft * x
        y_shaft = r_shaft * y

        z_transition =  1-head_fraction

        vertices_and_normals = []


        #
        # Head
        #

        rz = np.array([r_head, head_fraction])

        # Yes, they're swapped, because we want the normal perpendicular
        head_normal_z, head_normal_r = rz / np.sqrt(np.sum(rz**2))
        head_normal_z *= -1

        x_mid = head_normal_r * np.cos(mid_angles)
        y_mid = head_normal_r * np.sin(mid_angles)

        section = []
        for i in range(n_points-1):
            x_0 = x_head[i]
            x_1 = x_head[i+1]

            y_0 = y_head[i]
            y_1 = y_head[i+1]

            x_normal_0 = head_normal_r * x[i]
            x_normal_1 = head_normal_r * x[i+1]

            y_normal_0 = head_normal_r * y[i]
            y_normal_1 = head_normal_r * y[i+1]

            section += [
                [x_0, y_0, z_transition, x_normal_0, y_normal_0, head_normal_z],
                [x_1, y_1, z_transition, x_normal_1, y_normal_1, head_normal_z],
                [0,   0,   1,            x_mid[i],   y_mid[i],   head_normal_z],
            ]


        vertices_and_normals.append(section)

        #
        # Ring
        #

        section = []
        for i in range(n_points-1):
            x_inner_0 = x_shaft[i]
            x_inner_1 = x_shaft[i+1]

            y_inner_0 = y_shaft[i]
            y_inner_1 = y_shaft[i+1]

            x_outer_0 = x_head[i]
            x_outer_1 = x_head[i+1]

            y_outer_0 = y_head[i]
            y_outer_1 = y_head[i+1]

            section += [
                [x_inner_0, y_inner_0, z_transition, 0, 0, -1],
                [x_inner_1, y_inner_1, z_transition, 0, 0, -1],
                [x_outer_0, y_outer_0, z_transition, 0, 0, -1],
                [x_inner_1, y_inner_1, z_transition, 0, 0, -1],
                [x_outer_1, y_outer_1, z_transition, 0, 0, -1],
                [x_outer_0, y_outer_0, z_transition, 0, 0, -1],
            ]

        vertices_and_normals.append(section)

        #
        # Sides
        #

        section=[]
        for i in range(n_points-1):
            x_0 = x_shaft[i]
            x_1 = x_shaft[i+1]

            y_0 = y_shaft[i]
            y_1 = y_shaft[i+1]

            nx_0 = x[i]
            nx_1 = x[i+1]

            ny_0 = y[i]
            ny_1 = y[i+1]

            section += [
                [x_0, y_0, 0,     nx_0, ny_0, 0],
                [x_1, y_1, 0,     nx_1, ny_1, 0],
                [x_0, y_0, z_transition, nx_0, ny_0, 0],
                [x_1, y_1, 0,     nx_1, ny_1, 0],
                [x_1, y_1, z_transition, nx_1, ny_1, 0],
                [x_0, y_0, z_transition, nx_0, ny_0, 0],
            ]

        vertices_and_normals.append(section)

        #
        # Base
        #

        section = []
        for i in range(n_points-1):
            x_0 = x_shaft[i]
            x_1 = x_shaft[i + 1]

            y_0 = y_shaft[i]
            y_1 = y_shaft[i + 1]

            section += [
                [x_1, y_1, 0, 0, 0, -1],
                [x_0, y_0, 0, 0, 0, -1],
                [0,   0,   0, 0, 0, -1],
            ]

        vertices_and_normals.append(section)

        vertices_and_normals = [np.array(v_n, dtype=np.float32) for v_n in vertices_and_normals]

        for v_n in vertices_and_normals:
            self.add_vertex_normal_data(v_n)
