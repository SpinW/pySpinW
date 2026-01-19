import logging
from dataclasses import dataclass
from typing import Any

import numpy as np
from OpenGL.GL import *

logger = logging.Logger("pyspinw.gui.rendering.model")


class ModelChunk:
    def __enter__(self):
        pass

    def __exit__(self, exc_type, exc, traceback):
        pass


class Model:

    def __init__(self):
        self._vaos = []
        self._vbos = []
        self.lengths = []

    def add_vertex_normal_data(self, vertices_and_normals: np.ndarray):
        print("Creating vaos/vbos")

        try:

            vao = glGenVertexArrays(1)
            vbo = glGenBuffers(1)

            glBindVertexArray(vao)
            glBindBuffer(GL_ARRAY_BUFFER, vbo)
            glBufferData(GL_ARRAY_BUFFER, vertices_and_normals.nbytes, vertices_and_normals, GL_STATIC_DRAW)

            stride = 6 * 4  # 6 floats, 4 bytes each

            # position -> location 0
            glVertexAttribPointer(
                0, 3, GL_FLOAT, GL_FALSE,
                stride, ctypes.c_void_p(0)
            )
            glEnableVertexAttribArray(0)

            # normal -> location 1
            glVertexAttribPointer(
                1, 3, GL_FLOAT, GL_FALSE,
                stride, ctypes.c_void_p(3 * 4)
            )
            glEnableVertexAttribArray(1)

            self._vaos.append(vao)
            self._vbos.append(vbo)
            self.lengths.append(len(vertices_and_normals))


            glBindBuffer(GL_ARRAY_BUFFER, 0)
            glBindVertexArray(0)

        except Exception as e:
            logger.error(e)

        print("done creating vaos/vbos")

    @property
    def vaos(self):
        return self._vaos

    @property
    def vbos(self):
        return self._vbos

    def render_triangles(self):
        for vao, length in zip(self.vaos, self.lengths):
            glBindVertexArray(vao)
            glDrawArrays(GL_TRIANGLES, 0, length)
            glBindVertexArray(0)