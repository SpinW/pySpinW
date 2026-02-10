import numpy as np

from OpenGL.GL import *


class WireframeCube:
    def __init__(self):
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

        self.vao = glGenVertexArrays(1)
        self.vbo = glGenBuffers(1)

        glBindVertexArray(self.vao)

        glBindBuffer(GL_ARRAY_BUFFER, self.vbo)
        glBufferData(GL_ARRAY_BUFFER, edges.nbytes, edges, GL_STATIC_DRAW)

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 12, ctypes.c_void_p(0)) # 3 floats @ 4 bytes = 12
        glEnableVertexAttribArray(0)

        glBindVertexArray(0)

    def render_wireframe(self):
        glBindVertexArray(self.vao)
        glLineWidth(1.0)
        glDrawArrays(GL_LINES, 0, 24) # 24 vertices
