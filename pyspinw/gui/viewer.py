""" Viewer window """
import sys

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QSplitter, QWidget, QVBoxLayout, QTextEdit, QApplication

from pyspinw.gui.crystalview import CrystalViewerWidget
from pyspinw.gui.rendermodel import RenderModel
from pyspinw.gui.renderoptions import DisplayOptions, DisplayOptionsToolbar
from pyspinw.gui.textdisplay import TextDisplay
from pyspinw.hamiltonian import Hamiltonian


class Viewer(QWidget):
    """ Main viewer class """

    def __init__(self, hamiltonian: Hamiltonian, parent=None):

        super().__init__(parent)

        render_model = RenderModel(hamiltonian)

        layout = QVBoxLayout()

        self.toolbar = DisplayOptionsToolbar()

        splitter = QSplitter(Qt.Horizontal)
        self.viewer = CrystalViewerWidget(render_model)
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

        #
        # Wire up text and graphics
        #

        self.text_display.hoverChanged.connect(self.on_text_hover_changed)
        self.text_display.selectionChanged.connect(self.on_text_selection_changed)
        self.viewer.hoverChanged.connect(self.on_render_hover_changed)
        self.viewer.selectionChanged.connect(self.on_render_selection_changed)


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

        super().closeEvent(event)


def show_hamiltonian(hamiltonian):
    """ Show a Hamiltonian in the viewer"""
    app = QApplication()

    # Useful for checking particular display options
    # app.styleHints().setColorScheme(Qt.ColorScheme.Dark)
    # app.styleHints().setColorScheme(Qt.ColorScheme.Light)

    widget = Viewer(hamiltonian)
    widget.resize(800, 600)
    widget.show()
    app.exec()
