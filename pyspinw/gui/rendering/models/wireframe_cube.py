"""Geometry for a wireframe cube."""

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


class WireframeCube:
    """Geometry data for a wireframe cube.

    Attributes
    ----------
    vao : int
        OpenGL vertex array object handle.
    vbo : int
        OpenGL vertex buffer object handle.
    """

    def __init__(self, webgl=None):
        edges = np.array([
            [0,0,0], [0,0,1], # Z direction edges
            [1,0,0], [1,0,1],
            [0,1,0], [0,1,1],
            [1,1,0], [1,1,1],
            [0,0,0], [0,1,0], # Y direction edges
            [1,0,0], [1,1,0],
            [0,0,1], [0,1,1],
            [1,0,1], [1,1,1],
            [0,0,0], [1,0,0], # X direction edges
            [0,1,0], [1,1,0],
            [0,0,1], [1,0,1],
            [0,1,1], [1,1,1],
        ], dtype=np.float32)
        if webgl is not None and HASWEBGL:
            edges = np.ascontiguousarray(edges.reshape(-1, 3), dtype=np.float32)
            self._webgl = webgl[0]
            vbo = webgl[0].create_buffer_ext(target="ARRAY_BUFFER", src_data=edges, usage="STATIC_DRAW",
                                             auto_execute=False)
            self.vao = webgl[0].create_vertex_array_ext(webgl[1], [(vbo, "3f32", 0)], auto_execute=False)
            self.vertex_count = int(edges.shape[0])
            return
        elif not HASOPENGL:
            raise RuntimeError('OpenGL library cannot be loaded and WebGL not specified')
        else:
            self._webgl = None

        self.vao = glGenVertexArrays(1)
        self.vbo = glGenBuffers(1)

        glBindVertexArray(self.vao)

        glBindBuffer(GL_ARRAY_BUFFER, self.vbo)
        glBufferData(GL_ARRAY_BUFFER, edges.nbytes, edges, GL_STATIC_DRAW)

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 12, ctypes.c_void_p(0)) # 3 floats @ 4 bytes = 12
        glEnableVertexAttribArray(0)

        glBindVertexArray(0)

    def render_wireframe(self):
        """Render the wireframe using GL_LINES."""
        if self._webgl is not None:
            self._webgl.bind_vertex_array(self.vao)
            self._webgl.draw_arrays("LINES", 0, self.vertex_count)
            self._webgl.bind_vertex_array(None)
            return
        glBindVertexArray(self.vao)
        glLineWidth(1.0)
        glDrawArrays(GL_LINES, 0, 24) # 24 vertices
