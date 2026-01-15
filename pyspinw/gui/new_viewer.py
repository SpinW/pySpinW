import numpy as np
from PySide6.QtCore import QTimer
from PySide6.QtOpenGLWidgets import QOpenGLWidget
from PySide6.QtWidgets import QApplication
from OpenGL.GL import *
import sys

from pyspinw.gui.camera import Camera
from pyspinw.gui.load_shaders import load_shaders


class HamiltonianRenderer:
    def set_hamiltonian(self):
        pass

    def initialize(self):
        pass

    def render(self):
        pass

# Cube vertices
vertices = np.array([
    # front
    -1, -1,  1,  1, -1,  1,  1,  1,  1,
    -1, -1,  1,  1,  1,  1, -1,  1,  1,
    # back
    -1, -1, -1, -1,  1, -1,  1,  1, -1,
    -1, -1, -1,  1,  1, -1,  1, -1, -1,
    # left
    -1, -1, -1, -1, -1,  1, -1,  1,  1,
    -1, -1, -1, -1,  1,  1, -1,  1, -1,
    # right
     1, -1, -1,  1,  1, -1,  1,  1,  1,
     1, -1, -1,  1,  1,  1,  1, -1,  1,
    # top
    -1,  1, -1, -1,  1,  1,  1,  1,  1,
    -1,  1, -1,  1,  1,  1,  1,  1, -1,
    # bottom
    -1, -1, -1,  1, -1, -1,  1, -1,  1,
    -1, -1, -1,  1, -1,  1, -1, -1,  1,
], dtype=np.float32)


class GLWidget(QOpenGLWidget):
    def __init__(self):

        super().__init__()

        self.angle = 0.0
        self.radius = 10.0
        self.camera = Camera()
        self.shader_program = load_shaders()

    def initializeGL(self):
        glEnable(GL_DEPTH_TEST)

        self.vao = glGenVertexArrays(1)
        self.vbo = glGenBuffers(1)

        glBindVertexArray(self.vao)
        glBindBuffer(GL_ARRAY_BUFFER, self.vbo)
        glBufferData(GL_ARRAY_BUFFER, vertices.nbytes, vertices, GL_STATIC_DRAW)

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, None)
        glEnableVertexAttribArray(0)

        glBindBuffer(GL_ARRAY_BUFFER, 0)
        glBindVertexArray(0)

        self.angle = 0.0

        # Animation timer
        self.timer = QTimer()
        self.timer.timeout.connect(self.update)
        self.timer.start(16)


    def paintGL(self):
        glClearColor(0.05, 0.05, 0.08, 1.0)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        self.angle += 0.01

        camera_x = self.radius * np.sin(self.angle)
        camera_y = self.radius * np.cos(self.angle)


        self.camera.position = (camera_x, camera_y, 4)


        view = self.camera.view_matrix()
        proj = self.camera.perspective_matrix(0.01, 100)
        mvp = proj @ view

        glUseProgram(self.shader_program)
        loc = glGetUniformLocation(self.shader_program, "mvp")
        glUniformMatrix4fv(loc, 1, GL_FALSE, mvp.T)

        glBindVertexArray(self.vao)
        glDrawArrays(GL_TRIANGLES, 0, 36)
        glBindVertexArray(0)

    def resizeGL(self, w, h):
        self.camera.horizontal_pixels = w
        self.camera.vertical_pixels = h
        glViewport(0, 0, w, h)

    def mousePressEvent(self, event):
        pos = event.position()

    def mouseMoveEvent(self, event):
        pass

app = QApplication(sys.argv)
widget = GLWidget()
widget.resize(800, 600)
widget.show()
sys.exit(app.exec())

