""" Main widget for Qt/GL rendering of magnetic crystal structures"""

import numpy as np
from PySide6.QtCore import QTimer, QPoint
from PySide6.QtGui import Qt, QKeyEvent
from PySide6.QtOpenGLWidgets import QOpenGLWidget
from PySide6.QtWidgets import QApplication
from OpenGL.GL import *
import sys

from pyspinw.gui.buffers import IntegerBuffer
from pyspinw.gui.camera import Camera
from pyspinw.gui.render_model import RenderModel
from pyspinw.gui.rendering.models.arrow import Arrow
from pyspinw.gui.rendering.models.sphere import Sphere
from pyspinw.gui.rendering.models.tube import Tube

import logging

from pyspinw.gui.rendering.models.wrireframe_cube import WireframeCube
from pyspinw.gui.rendering.shader import SelectionShader, ObjectShader, CellShader, IDShader
from pyspinw.gui.renderoptions import DisplayOptions
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

        self.display_options = DisplayOptions()

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

        # These variables are used for selection / highlighting
        self.mouse_position: QPoint | None = None
        self.hover_index: int = 0

        # These variables are used for handling mouse dragging
        self.mouse_data: tuple[QPoint, Qt.MouseButton] | None= None
        self.mouse_rotation = np.eye(3)
        self.mouse_origin = np.array([0,0,0], dtype=float)

        # We want to be able to do hovering, and respond to key presses
        self.setMouseTracking(True)
        self.setFocusPolicy(Qt.StrongFocus)


    def initializeGL(self):
        """ Qt override, set up the GL rendering """
        glEnable(GL_DEPTH_TEST)

        try:
            # Normal objects
            self.sphere= Sphere(3)
            self.tube = Tube()
            self.arrow = Arrow()
            self.cube = WireframeCube()

            # Normal shaders
            self.object_shader = ObjectShader()
            self.selection_shader = SelectionShader()
            self.cell_shader = CellShader()

            self.id_shader = IDShader()

            # Framebuffer for ID rendering
            self.id_framebuffer = IntegerBuffer()


            self.angle = 0.0

            # Animation timer
            self.timer = QTimer()
            self.timer.timeout.connect(self.update)
            self.timer.start(16)


        except Exception as e:
            logger.exception(e)

    def paintGL(self):
        """ Qt override, paints the GL canvas"""

        #
        # Normal painting, let Qt do its thing
        #

        glEnable(GL_DEPTH_TEST)
        glEnable(GL_CULL_FACE)
        glDisable(GL_BLEND)

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

        moment_scale = 2 * self.display_options.atom_moment_scaling
        moment_scale_matrix = np.diag([moment_scale, moment_scale, moment_scale, 1])

        coupling_scale = 0.1 * self.display_options.coupling_scaling
        coupling_scaling = np.diag([coupling_scale, coupling_scale, 1, 1])

        if self.object_shader is not None and self.selection_shader is not None:

            self.object_shader.camera = self.camera
            self.selection_shader.camera = self.camera

            # Sites
            if self.display_options.show_sites:
                self.object_shader.object_color = 0.7, 0.8, 0.6


                for site in self.render_model.sites:

                    scaled_matrix = site.model_matrix @ moment_scale_matrix

                    if site.render_id == self.hover_index:
                        self.selection_shader.model_matrix = scaled_matrix
                        self.selection_shader.use()
                        self.arrow.render_back_wireframe()

                    self.object_shader.model_matrix = scaled_matrix
                    self.object_shader.use()
                    self.arrow.render_triangles()

            # Couplings
            if self.display_options.show_couplings:

                self.object_shader.object_color = 0.2, 0.4, 0.8

                for coupling in self.render_model.couplings:

                    scaled_matrix = coupling.model_matrix @ coupling_scaling

                    if coupling.render_id == self.hover_index:
                        self.selection_shader.model_matrix = scaled_matrix
                        self.selection_shader.use()
                        self.arrow.render_back_wireframe()

                    self.object_shader.model_matrix = scaled_matrix
                    self.object_shader.use()
                    self.tube.render_triangles()



        # Draw the supercell
        if self.display_options.show_supercell:
            supercell_matrix = np.eye(4, dtype=np.float32)
            supercell_matrix[:3, :3] = self.render_model.expanded.structure.unit_cell._xyz.T


            self.cell_shader.model_matrix = supercell_matrix
            self.cell_shader.camera = self.camera

            self.cell_shader.use()
            self.cube.render_wireframe()

        # Render unit cells
        if self.display_options.show_unit_cell:
            unit_cell_transform = np.eye(4, dtype=np.float32)
            unit_cell_transform[:3,:3] = self.render_model.hamiltonian.structure.unit_cell._xyz.T

            self.cell_shader.camera = self.camera

            for offset in self.render_model.hamiltonian.structure.supercell.cells():
                translation_matrix = np.eye(4, dtype=np.float32)
                translation_matrix[:3, 3] = offset.vector

                self.cell_shader.model_matrix = unit_cell_transform @ translation_matrix

                self.cell_shader.use()
                self.cube.render_wireframe()

        #
        # ID framebuffer
        #

        self.id_framebuffer.use(self.width(), self.height())

        if self.id_shader is not None:

            self.id_shader.camera = self.camera

            # Sites
            if self.display_options.show_sites:

                for site in self.render_model.sites:

                    scaled_matrix = site.model_matrix @ moment_scale_matrix

                    self.id_shader.model_matrix = scaled_matrix
                    self.id_shader.id_value = site.render_id
                    self.id_shader.use()
                    self.arrow.render_triangles()

            # Couplings
            if self.display_options.show_couplings:

                for coupling in self.render_model.couplings:
                    self.id_shader.id_value = coupling.render_id
                    self.id_shader.model_matrix = coupling.model_matrix @ coupling_scaling
                    self.id_shader.use()
                    self.tube.render_triangles()

        #
        # Calculate the object for the mouse over, and set its hover state
        #

        if self.mouse_position is not None:

            id = np.zeros(1, dtype=np.uint32) # Buffer to set data in

            x, y = self.mouse_position.x(), self.height() - self.mouse_position.y()

            glReadPixels(
                x, y,
                1, 1,
                GL_RED_INTEGER,
                GL_UNSIGNED_INT,
                id
            )

            # print(f"Hovering over ({x}, {y}) ID={id}")

            self.hover_index = int(id)

        else:
            self.hover_index = 0



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

    def leaveEvent(self, event):
        """ Qt override, mouse leaves widget"""
        self.mouse_position = None
        event.accept()

    def mouseMoveEvent(self, event):
        """Qt override, called on mouse movement"""


        self.mouse_position = event.position()

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


    def keyPressEvent(self, event: QKeyEvent):
        # We have this for debugging

        if event.key() == Qt.Key_W:

            import matplotlib.pyplot as plt
            im = self.id_framebuffer.get_image(self.width(), self.height())

            plt.figure()
            plt.imshow(im)
            plt.show()




        elif event.key() == Qt.Key_Escape:
            print("Escape pressed")

        event.accept()