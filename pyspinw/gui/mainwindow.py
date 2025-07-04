""" Entrypoint and main window class for the GUI"""

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QApplication

from pyspinw.gui.crystalviewer.viewer import CrystalViewer
from pyspinw.gui.symmetryeditor import SymmetryEditor
from pyspinw.gui.sitecouplingeditor import SiteAndCouplingEditor
from pyspinw.symmetry.unitcell import UnitCell


class MainWindow(QMainWindow):
    """ Main window """

    def __init__(self):
        super().__init__()

        self.setWindowTitle("pySpinW Crystal Editor")

        self.viewer = CrystalViewer(parent=self)
        self.lattice_parameters = SymmetryEditor()
        self.sites = SiteAndCouplingEditor()

        self.setCentralWidget(self.viewer)

        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self.lattice_parameters)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.sites)

        # Connections
        self.lattice_parameters.unit_cell_widget.unit_cell_changed.connect(self._on_unit_cell_changed)
        self.lattice_parameters.symmetry_changed.connect(self._on_symmetry_changed)
        self.sites.graphics_relevant_change.connect(self._on_graphics_relevant_site_change)

    def _on_unit_cell_changed(self, unit_cell: UnitCell):
        self.viewer.unit_cell = unit_cell

    def _on_symmetry_changed(self):
        self.sites.symmetry = self.lattice_parameters.symmetry

    def _on_graphics_relevant_site_change(self):
        self.viewer.sites = self.sites.sites_for_drawing
        self.viewer.update()

if __name__ == "__main__":
    app = QApplication([])

    editor = MainWindow()
    editor.viewer.unit_cell = UnitCell(1,1,1,90,90,90)
    editor.show()

    app.exec_()
