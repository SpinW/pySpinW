""" Main widget for Qt/GL rendering of magnetic crystal structures"""

import numpy as np
from PySide6.QtCore import QTimer, QPoint
from PySide6.QtGui import Qt
from PySide6.QtOpenGLWidgets import QOpenGLWidget
from PySide6.QtWidgets import QApplication
from OpenGL.GL import *
import sys

from pyspinw.gui.camera import Camera
from pyspinw.gui.render_model import RenderModel
from pyspinw.gui.rendering.models.arrow import Arrow
from pyspinw.gui.rendering.models.sphere import Sphere
from pyspinw.gui.rendering.models.tube import Tube

import logging

from pyspinw.gui.rendering.models.wrireframe_cube import WireframeCube
from pyspinw.gui.rendering.shader import SelectionShader, ObjectShader, CellShader
from pyspinw.util import rotation_matrix

logger = logging.Logger(__name__)



class CrystalViewerWidget(QOpenGLWidget):
    """ Qt widget to show magnetic crystal structures """

    mouse_angle_sensitivity = 0.01
    mouse_move_sensitivity = 0.005
    mouse_zoom_sensitivity = 0.002
    min_view_radius = 0.01

    def __init__(self, render_model: RenderModel):

        super().__init__()

        self.render_model = render_model

        self.camera = Camera()
        # self.shader_program = None
        # self.default_shader = None
        #
        self.object_shader: ObjectShader | None = None
        self.selection_shader: SelectionShader | None = None

        # Set up antialiasing
        format = self.format()
        format.setSamples(4)
        self.setFormat(format)

        # View details
        self.view_origin = render_model.expanded.structure.unit_cell.centre #np.array([0,0,0], dtype=float)
        self.view_rotation = np.eye(3)
        self.view_radius = 10.0


        # These variables are used for handling mouse dragging
        self.mouse_data: tuple[QPoint, Qt.MouseButton] | None= None
        self.mouse_rotation = np.eye(3)
        self.mouse_origin = np.array([0,0,0], dtype=float)



    def initializeGL(self):
        """ Qt override, set up the GL rendering """
        glEnable(GL_DEPTH_TEST)

        try:

            self.sphere1 = Sphere(3)
            self.sphere2 = Sphere(3)
            self.tube = Tube()
            self.arrow = Arrow()
            self.cube = WireframeCube()

            # self.shader_program = load_shaders(vertex_filename="phong_vertex", fragment_filename="tailored_fragment")
            # self.default_shader = load_shaders()
            # self.shader_program = load_shaders(vertex_filename="phong_vertex", fragment_filename="default_fragment")
            # self.shader_program = load_shaders()

            self.object_shader = ObjectShader()
            self.selection_shader = SelectionShader()
            self.cell_shader = CellShader()

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

        # Work out camera movement stuff

        compound_rotation = self.view_rotation @ self.mouse_rotation

        camera_world = self.view_origin + self.mouse_origin + \
                       compound_rotation @ np.array([0, 0, self.view_radius], dtype=float)

        origin_world = self.view_origin + self.mouse_origin

        up_world = compound_rotation @ np.array([0,1,0], dtype=float)

        self.camera.position = tuple(camera_world)
        self.camera.up = tuple(up_world)
        self.camera.look_at = tuple(origin_world)

        # Do the rendering

        if self.object_shader is not None and self.selection_shader is not None:
            self.object_shader.object_color = 0.7, 0.8, 0.6

            self.object_shader.camera = self.camera
            self.selection_shader.camera = self.camera


            for site in self.render_model.sites:

                if site.is_selected:
                    self.selection_shader.model_matrix = site.model_matrix
                    self.selection_shader.use()
                    self.arrow.render_back_wireframe()

                self.object_shader.model_matrix = site.model_matrix
                self.object_shader.use()
                self.arrow.render_triangles()

            self.object_shader.object_color = 0.2, 0.4, 0.8
            thickness = 0.02
            coupling_scaling = np.diag([thickness, thickness, 1, 1])
            for coupling in self.render_model.couplings:
                self.object_shader.model_matrix = coupling.model_matrix @ coupling_scaling
                self.object_shader.use()
                self.tube.render_triangles()

        # Draw the supercell
        supercell_matrix = np.eye(4, dtype=np.float32)
        supercell_matrix[:3, :3] = self.render_model.expanded.structure.unit_cell._xyz.T


        self.cell_shader.model_matrix = supercell_matrix
        self.cell_shader.camera = self.camera

        self.cell_shader.use()
        self.cube.render_wireframe()

        # Render unit cells

        unit_cell_transform = np.eye(4, dtype=np.float32)
        unit_cell_transform[:3,:3] = self.render_model.hamiltonian.structure.unit_cell._xyz.T
        for offset in self.render_model.hamiltonian.structure.supercell.cells():
            translation_matrix = np.eye(4, dtype=np.float32)
            translation_matrix[:3, 3] = offset.vector

            self.cell_shader.model_matrix = unit_cell_transform @ translation_matrix

            self.cell_shader.use()
            self.cube.render_wireframe()



    def resizeGL(self, w, h):
        """ Qt override, called when window is resized """
        self.camera.horizontal_pixels = w
        self.camera.vertical_pixels = h
        glViewport(0, 0, w, h)

    def mousePressEvent(self, event):
        """ Qt override, called on mouse press"""

        if self.mouse_data is None:
            self.mouse_data = event.position(), event.button()

            if event.button() == Qt.MouseButton.LeftButton:
                self.mouse_rotation = np.eye(3)

            elif event.button() == Qt.MouseButton.RightButton:
                self.mouse_rotation = np.eye(3)
                self.mouse_origin = np.zeros((3, ), dtype=float)

    def reset_view(self):
        """ Reset the view to default """
        self.view_rotation = np.eye(3)
        self.view_origin = np.zeros((3, ), dtype=float)

    def mouseReleaseEvent(self, event):
        """ Qt override, called on mouse up"""

        self.mouse_data = None

        # self.view_rotation = self.mouse_rotation @ self.view_rotation
        self.view_rotation =  self.view_rotation @ self.mouse_rotation
        self.view_origin = self.mouse_origin + self.view_origin

        self.mouse_rotation = np.eye(3, dtype=float)
        self.mouse_origin = np.zeros((3, ), dtype=float)

    def wheelEvent(self, event):
        """ Qt override, called when scroll wheel moves """
        self.view_radius *= 2**(self.mouse_zoom_sensitivity * event.angleDelta().y())

        if self.view_radius < self.min_view_radius:
            self.view_radius = self.min_view_radius

    def mouseMoveEvent(self, event):
        """Qt override, called on mouse movement"""

        if self.mouse_data is not None:
            start, button = self.mouse_data
            end = event.position()

            dx = end.x() - start.x()
            dy = end.y() - start.y()

            if button == Qt.MouseButton.LeftButton:
                self.mouse_rotation = \
                    rotation_matrix(-self.mouse_angle_sensitivity*dx, (0, 1, 0)) @ \
                     rotation_matrix(-self.mouse_angle_sensitivity*dy, (1, 0, 0))


            if button == Qt.MouseButton.RightButton:
                #r = np.linalg.inv(self.view_rotation)
                r = self.view_rotation
                self.mouse_origin = r @ np.array([
                    -self.mouse_move_sensitivity*dx,
                    self.mouse_move_sensitivity*dy,
                    0.0])

def show_hamiltonian(render_model):



    app = QApplication()
    widget = CrystalViewerWidget(render_model)
    widget.resize(800, 600)
    widget.show()
    sys.exit(app.exec())

