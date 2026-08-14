""" Camera model """

import numpy as np


class Camera:
    """ Data defining a camera """

    def __init__(self,
        position: tuple[float, float, float] = (-10, 0, 0),
        look_at: tuple[float, float, float] = (0,0,0),
        up: tuple[float, float, float] = (0,0,1),
        fov_deg: float = 50.0,
        horizontal_pixels: int = 800,
        vertical_pixels: int = 600,
        orthographic=True):

        self.position = position
        self.look_at = look_at
        self.up = up
        self.fov_deg = fov_deg
        self.horizontal_pixels = horizontal_pixels
        self.vertical_pixels = vertical_pixels
        self.orthographic = orthographic
        self.fixed_size: float | None = None

    @property
    def aspect_ratio(self) -> float:
        """ Aspect ratio of this camera """
        return self.horizontal_pixels / self.vertical_pixels

    @property
    def distance(self) -> float:
        """ Distance from camera to view origin """
        return float(np.sqrt(np.sum((np.array(self.position) - np.array(self.look_at))**2)))

    def projection_matrix(self, near: float = 0.01, far: float = 100):
        """ Get the projection matrix (orthographic/perspective based on options)

        :param near: near clipping plane distance
        :param far: far clipping plane distance

        :return: 4x4 float32 frustum matrix
        """
        if self.orthographic:
            return self.orthographic_matrix(near, far)
        else:
            return self.perspective_matrix(near, far)

    def perspective_matrix(self, near: float = 0.01, far: float = 100):
        """ Get the perspective projection matrix for GL

        :param near: near clipping plane distance
        :param far: far clipping plane distance

        :return: 4x4 float32 frustum matrix
        """
        f = 1.0 / np.tan(np.radians(self.fov_deg) / 2)
        m = np.zeros((4, 4), dtype=np.float32)

        m[0, 0] = f / self.aspect_ratio
        m[1, 1] = f
        m[2, 2] = (far + near) / (near - far)
        m[2, 3] = (2 * far * near) / (near - far)
        m[3, 2] = -1

        return m

    def orthographic_matrix(self, near: float = 0.01, far: float = 100):
        """ Get the orthographic projection matrix for GL

        :param near: near clipping plane distance
        :param far: far clipping plane distance

        :return: 4x4 float32 frustum matrix
        """
        if self.fixed_size is None:
            height = 2.0 * self.distance * np.tan(np.radians(self.fov_deg) / 2)
        else:
            height = self.fixed_size

        width = height * self.aspect_ratio

        return np.array([
            [2/width, 0, 0, 0],
            [0, 2/height, 0, 0],
            [0, 0, 2/(near-far), (near+far)/(near-far)],
            [0, 0, 0, 1]
        ], dtype=np.float32)

    def view_matrix(self):
        """Transform from world position to camera relative """
        eye = np.array(self.position, dtype=np.float32)
        target = np.array(self.look_at, dtype=np.float32)
        up = np.array(self.up, dtype=np.float32)

        # Look at matrix
        f = target - eye
        f = f / np.linalg.norm(f)

        s = np.cross(f, up)
        s = s / np.linalg.norm(s)

        u = np.cross(s, f)

        view = np.identity(4, dtype=np.float32)
        view[0, :3] = s
        view[1, :3] = u
        view[2, :3] = -f
        view[0, 3] = -np.dot(s, eye)
        view[1, 3] = -np.dot(u, eye)
        view[2, 3] = np.dot(f, eye)

        return view


    axes_target = np.ones((3, ), dtype=np.float32) / 2
    def axes_view_matrix(self):
        """View matrix for axes"""
        eye = np.array(self.position, dtype=np.float32)
        eye_mag = np.sqrt(np.sum(eye ** 2))
        eye *= 3 / eye_mag
        up = np.array(self.up, dtype=np.float32)

        # Look at matrix
        f = self.axes_target - eye
        f = f / np.linalg.norm(f)

        s = np.cross(f, up)
        s = s / np.linalg.norm(s)

        u = np.cross(s, f)

        view = np.identity(4, dtype=np.float32)
        view[0, :3] = s
        view[1, :3] = u
        view[2, :3] = -f
        view[0, 3] = -np.dot(s, eye)
        view[1, 3] = -np.dot(u, eye)
        view[2, 3] = np.dot(f, eye)

        return view
