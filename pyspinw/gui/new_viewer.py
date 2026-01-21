""" Main widget for Qt/GL rendering of magnetic crystal structures"""

import numpy as np
from PySide6.QtCore import QTimer, QPoint
from PySide6.QtGui import QSurfaceFormat, Qt
from PySide6.QtOpenGLWidgets import QOpenGLWidget
from PySide6.QtWidgets import QApplication
from OpenGL.GL import *
import sys

from pyspinw.gui.camera import Camera
from pyspinw.gui.load_shaders import load_shaders
from pyspinw.gui.rendering.arrow import Arrow
from pyspinw.gui.rendering.sphere import Sphere
from pyspinw.gui.rendering.tube import Tube

import logging

logger = logging.Logger(__name__)



class CrystalViewerWidget(QOpenGLWidget):
    """ Qt widget to show magnetic crystal structures """

    def __init__(self):

        super().__init__()

        self.angle = 0.0
        self.radius = 10.0
        self.camera = Camera()
        self.shader_program = None

        # Set up antialiasing
        format = self.format()
        format.setSamples(4)
        self.setFormat(format)

        self.mouse_data: tuple[QPoint, Qt.MouseButton] | None= None

        self.view_origin = np.array([0,0,0], dtype=float)
        self.view_rotation = np.eye(3)



    def initializeGL(self):
        """ Qt override, set up the GL rendering """
        glEnable(GL_DEPTH_TEST)

        try:

            self.sphere1 = Sphere(3)
            self.sphere2 = Sphere(3)
            self.tube = Tube()
            self.arrow = Arrow()

            self.shader_program = load_shaders(vertex_filename="phong_vertex", fragment_filename="tailored_fragment")
            # self.shader_program = load_shaders(vertex_filename="phong_vertex", fragment_filename="default_fragment")
            # self.shader_program = load_shaders()

            self.angle = 0.0

            # Animation timer
            self.timer = QTimer()
            self.timer.timeout.connect(self.update)
            self.timer.start(16)

        except Exception as e:
            logger.exception(e)

    def paintGL(self):
        """ Qt override, paints the GL canvas"""
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

        # self.sphere1.render_triangles()
        # self.sphere2.render_triangles()

        # self.tube.render_triangles()
        self.arrow.render_triangles()

        #
        # glBindVertexArray(self.vao)
        # glDrawArrays(GL_TRIANGLES, 0, len(vertices_and_normals))
        # glBindVertexArray(0)

    def resizeGL(self, w, h):
        """ Qt override, called when window is resized """
        self.camera.horizontal_pixels = w
        self.camera.vertical_pixels = h
        glViewport(0, 0, w, h)

    def mousePressEvent(self, event):
        """ Qt override, called on mouse press"""

        if self.mouse_data is None:
            self.mouse_data = event.position(), event.button()


    def mouseReleaseEvent(self, event):
        """ Qt override, called on mouse up"""

        self.mouse_data = None


    def mouseMoveEvent(self, event):
        """Qt override, called on mouse movement"""

        if self.mouse_data is not None:


app = QApplication(sys.argv)
widget = CrystalViewerWidget()
widget.resize(800, 600)
widget.show()
sys.exit(app.exec())

