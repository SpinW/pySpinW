from PySide6.QtCore import Qt
from PySide6.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QApplication

from pyspinw.gui.crystalviewer.viewer import CrystalViewer
from pyspinw.gui.lattice import LatticeParameters
from pyspinw.unitcell import UnitCell


class CrystalEditor(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("pySpinW Crystal Editor")

        self.viewer = CrystalViewer(parent=self)
        self.lattice_parameters = LatticeParameters()

        self.setCentralWidget(self.viewer)

        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self.lattice_parameters)

        # Connections
        self.lattice_parameters.unit_cell_changed.connect(self._on_unit_cell_changed)

    def _on_unit_cell_changed(self, unit_cell: UnitCell):
        self.viewer.unit_cell = unit_cell


if __name__ == "__main__":
    app = QApplication([])

    editor = CrystalEditor()
    editor.viewer.unit_cell = UnitCell(1,1,1,90,90,90)
    editor.show()

    app.exec_()