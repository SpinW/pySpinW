import numpy as np
from PySide6.QtCore import QTimer
from PySide6.QtOpenGLWidgets import QOpenGLWidget
from PySide6.QtWidgets import QApplication
from OpenGL.GL import *
import sys

from pyspinw.gui.camera import Camera
from pyspinw.gui.load_shaders import load_shaders
from pyspinw.gui.rendering.sphere import Sphere


class HamiltonianRenderer:
    def set_hamiltonian(self):
        pass

    def initialize(self):
        pass

    def render(self):
        pass

vertices_and_normals = Sphere(3).vertices_and_normals

class GLWidget(QOpenGLWidget):
    def __init__(self):

        super().__init__()

        self.angle = 0.0
        self.radius = 10.0
        self.camera = Camera()
        self.shader_program = None

    def initializeGL(self):
        glEnable(GL_DEPTH_TEST)

        self.vao = glGenVertexArrays(1)
        self.vbo = glGenBuffers(1)

        glBindVertexArray(self.vao)
        glBindBuffer(GL_ARRAY_BUFFER, self.vbo)
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

        glBindBuffer(GL_ARRAY_BUFFER, 0)
        glBindVertexArray(0)

        self.shader_program = load_shaders(vertex_filename="phong_vertex", fragment_filename="phong_fragment")
        # self.shader_program = load_shaders(vertex_filename="phong_vertex", fragment_filename="default_fragment")
        # self.shader_program = load_shaders()

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

        camera_position = (camera_x, camera_y, 4)

        self.camera.position = camera_position


        view = self.camera.view_matrix()
        proj = self.camera.perspective_matrix(0.01, 100)
        projectionView = proj @ view

        glUseProgram(self.shader_program)

        model_loc = glGetUniformLocation(self.shader_program, "model")
        projection_view_loc = glGetUniformLocation(self.shader_program, "projectionView")

        light_pos_loc = glGetUniformLocation(self.shader_program, "lightPos")
        view_pos_loc = glGetUniformLocation(self.shader_program, "viewPos")
        light_color_loc = glGetUniformLocation(self.shader_program, "lightColor")
        object_color_loc = glGetUniformLocation(self.shader_program, "objectColor")

        glUniformMatrix4fv(model_loc, 1, GL_FALSE, np.eye(4, dtype=np.float32))
        glUniformMatrix4fv(projection_view_loc, 1, GL_FALSE, projectionView.T)

        glUniform3f(light_pos_loc, -20,0,0)
        glUniform3f(view_pos_loc, *camera_position)

        glUniform3f(light_color_loc, 1, 1, 1)
        glUniform3f(object_color_loc, 0.7, 0.8, 0.6)

        glBindVertexArray(self.vao)
        glDrawArrays(GL_TRIANGLES, 0, len(vertices_and_normals))
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

