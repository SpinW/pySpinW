"""Base class for renderable objects."""

import numpy as np

# Must have one of either ipywebgl or OpenGL
try:
    import ipywebgl
    HASWEBGL = True
except ModuleNotFoundError:
    HASWEBGL = False

try:
    from OpenGL.GL import *
    HASOPENGL = True
except ModuleNotFoundError:
    HASOPENGL = False
    if not HASWEBGL:
        raise

import logging
logger = logging.Logger(__name__)

class Model:
    """Base class for renderable objects."""

    def __init__(self, webgl=None):
        self._vaos = []
        self._vbos = []
        self._lengths = []
        self._webgldata = []
        if webgl is not None and HASWEBGL:
            self._webgl, self._webglprogram = webgl # tuple input
        elif not HASOPENGL:
            raise RuntimeError('OpenGL library cannot be loaded and WebGL not specified')
        else:
            self._webgl = None

    def add_vertex_normal_data(self, vertices_and_normals: np.ndarray):
        """Add vertex and normal data.

        Requires an n-by-6 numpy array of float32s.

        Parameters
        ----------
        vertices_and_normals : (n, 6) numpy.ndarray
            Array of interleaved 3D vertex positions and 3D normal vectors of
            dtype ``float32``.
        """
        if self._webgl is not None:
            self._webgldata.append(vertices_and_normals)
            return

        vao = glGenVertexArrays(1)
        vbo = glGenBuffers(1)

        glBindVertexArray(vao)
        glBindBuffer(GL_ARRAY_BUFFER, vbo)
        glBufferData(GL_ARRAY_BUFFER, vertices_and_normals.nbytes, vertices_and_normals, GL_STATIC_DRAW)

        stride = 6 * 4

        glVertexAttribPointer(
            0, 3, GL_FLOAT, GL_FALSE,
            stride, ctypes.c_void_p(0)
        )
        glEnableVertexAttribArray(0)

        glVertexAttribPointer(
            1, 3, GL_FLOAT, GL_FALSE,
            stride, ctypes.c_void_p(3 * 4)
        )
        glEnableVertexAttribArray(1)

        self._vaos.append(vao)
        self._vbos.append(vbo)
        self._lengths.append(len(vertices_and_normals))

        glBindBuffer(GL_ARRAY_BUFFER, 0)
        glBindVertexArray(0)


    def _render_webgldata(self, face):
        if isinstance(self._webgldata, list):
            self._webgldata = np.ascontiguousarray(np.vstack(self._webgldata), dtype=np.float32)
            self.vertex_count = int(self._webgldata.shape[0])
            vbo = self._webgl.create_buffer_ext(target="ARRAY_BUFFER", src_data=self._webgldata, usage="STATIC_DRAW",
                                                auto_execute=False)
            self.vao = self._webgl.create_vertex_array_ext(self._webglprogram, [(vbo, "3f32 3f32", 0, 1)],
                                                           auto_execute=False)
        self._webgl.enable(cull_face=True)
        self._webgl.cull_face(face)
        self._webgl.bind_vertex_array(self.vao)
        self._webgl.draw_arrays("TRIANGLES", 0, self.vertex_count)
        self._webgl.bind_vertex_array(None)


    def render_triangles(self):
        """Call the GL rendering for the triangles in this model."""
        if self._webgl is not None:
            return self._render_webgldata("BACK")
        glEnable(GL_CULL_FACE)
        glCullFace(GL_BACK)
        for vao, length in zip(self._vaos, self._lengths):
            glBindVertexArray(vao)
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
            glDrawArrays(GL_TRIANGLES, 0, length)
        glBindVertexArray(0)

    def render_back(self):
        """Render the back faces only - used for making selections."""
        if self._webgl is not None:
            return self._render_webgldata("FRONT")
        glEnable(GL_CULL_FACE)
        glCullFace(GL_FRONT)
        for vao, length in zip(self._vaos, self._lengths):
            glBindVertexArray(vao)
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
            glDrawArrays(GL_TRIANGLES, 0, length)


        glBindVertexArray(0)
