""" Base class for renderable objects """

import logging

import numpy as np
from OpenGL.GL import *

logger = logging.Logger(__name__)

class Model:
    """ Base class for renderable objects """

    def __init__(self):
        self._vaos = []
        self._vbos = []
        self._lengths = []

    def add_vertex_normal_data(self, vertices_and_normals: np.ndarray):
        """ Add vertex and normal data

        Requires an n-by-6 numpy array of float32s
        """
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


    def render_triangles(self):
        """ Call the GL rendering for the triangles in this model """
        glEnable(GL_CULL_FACE)
        glCullFace(GL_BACK)
        for vao, length in zip(self._vaos, self._lengths):
            glBindVertexArray(vao)
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
            glDrawArrays(GL_TRIANGLES, 0, length)
        glBindVertexArray(0)

    def render_back_wireframe(self):
        """ Render the back faces using wireframe - used for making selections """
        glEnable(GL_CULL_FACE)
        glCullFace(GL_FRONT)
        for vao, length in zip(self._vaos, self._lengths):
            glBindVertexArray(vao)
            glLineWidth(4)
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
            glDrawArrays(GL_TRIANGLES, 0, length)
        glBindVertexArray(0)
