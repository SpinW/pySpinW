""" Viewer window """
import sys

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QSplitter, QWidget, QVBoxLayout, QTextEdit, QApplication

from pyspinw.gui.new_crystalviewer import CrystalViewerWidget
from pyspinw.gui.render_model import RenderModel
from pyspinw.gui.renderoptions import DisplayOptions, DisplayOptionsToolbar
from pyspinw.gui.text_display import TextDisplay
from pyspinw.hamiltonian import Hamiltonian


class Viewer(QWidget):
    def __init__(self, hamiltonian: Hamiltonian, parent=None):

        super().__init__(parent)

        render_model = RenderModel(hamiltonian)

        splitter = QSplitter(Qt.Horizontal)

        self.viewer = CrystalViewerWidget(render_model)
        self.text_display = TextDisplay(render_model)

        splitter.addWidget(self.viewer)
        splitter.addWidget(self.text_display)

        layout = QVBoxLayout()

        self.toolbar = DisplayOptionsToolbar()

        layout.addWidget(self.toolbar)

        layout.addWidget(splitter, stretch=1)

        self.setLayout(layout)

        #
        # Wire up settings
        #

        self.viewer.display_options = self.toolbar.display_options()
        self.toolbar.displayOptionsChanged.connect(self.on_display_options_changed)

        #
        # Wire up left and right
        #

        self.text_display.hoverChanged.connect(self.on_text_hover_changed)

    def on_display_options_changed(self):
        self.viewer.display_options = self.toolbar.display_options()

    def on_text_hover_changed(self):
        self.viewer.hover_ids = self.text_display.hover_ids

    def closeEvent(self, event):
        self.toolbar.save_settings()
        super().closeEvent(event)


def show_hamiltonian(hamiltonian):

    app = QApplication()
    widget = Viewer(hamiltonian)
    widget.resize(800, 600)
    widget.show()
    sys.exit(app.exec())
