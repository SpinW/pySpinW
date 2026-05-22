""" Viewer window """
import ctypes
import os
import sys
import threading

import numpy as np
from imageio import imwrite

from PySide6.QtCore import Qt, QDir
from PySide6.QtWidgets import QSplitter, QWidget, QVBoxLayout, QTextEdit, QApplication, QFileDialog, QMessageBox
from numpy._typing import ArrayLike

from pyspinw import Structure
from pyspinw.gui.crystalview import CrystalViewerWidget
from pyspinw.gui.icons.iconload import png_icon
from pyspinw.gui.rendermodel import RenderModel
from pyspinw.gui.displayoptionstoolbar import DisplayOptionsToolbar
from pyspinw.gui.displayoptions import DisplayOptions
from pyspinw.gui.textdisplay import TextDisplay
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.util import rotation_matrix, rotation_from_z

_unique_id_counter = -1
def _generate_unique_id():
    """ Generate a unique ID for each site currently loaded"""
    global _unique_id_counter # noqa: PLW0603
    _unique_id_counter += 1
    return _unique_id_counter

class Viewer(QWidget):
    """ Main viewer class """

    def __init__(self,
                 hamiltonian: Hamiltonian,
                 initial_rotation: np.ndarray | None = None,
                 initial_distance: float | None = None,
                 display_options: DisplayOptions | None = None,
                 parent=None):

        super().__init__(parent)

        # Check parameters

        if initial_rotation is not None:
            if not np.allclose(initial_rotation @ initial_rotation.T, np.eye(3)):
                raise ValueError("Expected initial_rotation to be orthogonal")

            if not np.isclose(np.linalg.det(initial_rotation), 1):
                raise ValueError("Expected a rotation, not a rotoreflection")

        if initial_distance is not None and initial_distance <= 0:
            raise ValueError("Expected initial distance to be strictly positive")

        # Initialisation

        self._unique_id = _generate_unique_id()

        render_model = RenderModel(hamiltonian)

        layout = QVBoxLayout()

        self.toolbar = DisplayOptionsToolbar(initial_display_options=display_options)

        splitter = QSplitter(Qt.Horizontal)
        self.viewer = CrystalViewerWidget(render_model, initial_rotation, initial_distance)
        self.text_display = TextDisplay(render_model)
        splitter.addWidget(self.viewer)
        splitter.addWidget(self.text_display)

        layout.addWidget(self.toolbar)
        layout.addWidget(splitter, stretch=1)

        self.setLayout(layout)

        #
        # Wire up settings
        #

        self.viewer.display_options = self.toolbar.display_options()
        self.toolbar.displayOptionsChanged.connect(self.on_display_options_changed)
        self.toolbar.requestViewReset.connect(self.on_reset_view_requested)
        self.toolbar.requestSnapshot.connect(self.on_snapshot_requested)

        #
        # Wire up text and graphics
        #

        self.text_display.hoverChanged.connect(self.on_text_hover_changed)
        self.text_display.selectionChanged.connect(self.on_text_selection_changed)
        self.viewer.hoverChanged.connect(self.on_render_hover_changed)
        self.viewer.selectionChanged.connect(self.on_render_selection_changed)

    def on_snapshot_requested(self):
        """ Take a snapshot """
        dialog = QFileDialog(self)
        dialog.setAcceptMode(QFileDialog.AcceptSave)
        dialog.setOption(QFileDialog.DontConfirmOverwrite, True)
        dialog.setNameFilters([
            "PNG Images (*.png)",
            "JPEG Images (*.jpg *.jpeg)",
            "All Files (*)"
        ])
        dialog.selectNameFilter("PNG Images (*.png)")
        dialog.setDefaultSuffix("png")
        dialog.setDirectory(QDir.currentPath())


        if dialog.exec():
            filename = dialog.selectedFiles()[0]
        else:
            return

        if not filename:
            # Cancelled
            return

        if os.path.exists(filename):
            # File exists

            result = QMessageBox.question(
                self,
                "Overwrite file?",
                "This file already exists, do you want to replace it",
                QMessageBox.Ok | QMessageBox.Cancel
            )

            if not result == QMessageBox.Ok:
                return

        self.save_snapshot(filename)


    def save_snapshot(self, filename):
        """ Write the current viewer data to a file"""
        data = self.viewer.snapshot()
        imwrite(filename, data)


    def on_display_options_changed(self):
        """ Called when the toolbar changes"""
        self.viewer.display_options = self.toolbar.display_options()

    def on_text_hover_changed(self):
        """ Called when mouse hovers over text components """
        self.viewer.hover_ids = self.text_display.hover_ids

    def on_text_selection_changed(self):
        """ Called when text components are selected """
        self.viewer.current_selection = self.text_display.current_selection

    def on_render_hover_changed(self):
        """ Called when hovering over graphical components"""
        self.text_display.set_hover(self.viewer.hover_ids)

    def on_render_selection_changed(self):
        """ Called when selecting through the graphical window"""
        self.text_display.set_selection(self.viewer.current_selection)

    def on_reset_view_requested(self):
        """ Called when a view reset is requested"""
        self.viewer.reset_view()

    def closeEvent(self, event):
        """ Qt override for close event"""
        # Save settings on exit
        self.toolbar.save_settings()

        try:
            del _VIEWERS[self._unique_id]
        except Exception:
            pass

        super().closeEvent(event)


_APP = None
_VIEWERS = {}

def get_app():
    """ Get the Qt instance """
    global _APP #noqa: PLW0603
    app = QApplication.instance()
    if app is None:
        _APP = QApplication(sys.argv)
        app = _APP
    return app


def snapshot(object: Hamiltonian | Structure,
             view_point: ArrayLike = (0,0,10.0),
             display_options: DisplayOptions | None = None):

    """ Make (and display) an image of a structure/hamiltonian """

    try:
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("org.spinw.pyspinw")
    except Exception:
        pass

    if isinstance(object, Structure):
        object = Hamiltonian(object, [])

    if not isinstance(object, Hamiltonian):
        raise TypeError("Viewer needs to be given a Hamiltonian or Structure")

    # Parameters for setting the view
    view_point = np.array(view_point)
    rotation = rotation_from_z(view_point)
    distance = np.sqrt(np.sum(view_point**2))


    app = get_app()

    app.setWindowIcon(png_icon("pyspinw"))

    # Useful for checking particular display options
    # app.styleHints().setColorScheme(Qt.ColorScheme.Dark)
    # app.styleHints().setColorScheme(Qt.ColorScheme.Light)

    viewer = Viewer(object, rotation, distance, display_options)
    viewer.setWindowTitle("Hamiltonian Viewer")
    viewer.resize(800, 600)
    viewer.show()

    # Save a reference
    _VIEWERS[viewer._unique_id] = viewer


    app.exec()


def show_object(object: Hamiltonian | Structure, block=True):
    """ Show a Hamiltonian or structure in the viewer"""
    try:
        ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("org.spinw.pyspinw")
    except Exception:
        pass

    if isinstance(object, Structure):
        object = Hamiltonian(object, [])

    if not isinstance(object, Hamiltonian):
        raise TypeError("Viewer needs to be given a Hamiltonian or Structure")

    app = get_app()

    app.setWindowIcon(png_icon("pyspinw"))

    # Useful for checking particular display options
    # app.styleHints().setColorScheme(Qt.ColorScheme.Dark)
    # app.styleHints().setColorScheme(Qt.ColorScheme.Light)

    viewer = Viewer(object)
    viewer.setWindowTitle("Hamiltonian Viewer")
    viewer.resize(800, 600)
    viewer.show()

    # Save a reference
    _VIEWERS[viewer._unique_id] = viewer

    if block:
        app.exec()
