from PySide6.QtCore import Signal
from PySide6.QtWidgets import QBoxLayout, QGridLayout, QWidget, QApplication, QLabel, QDockWidget, QVBoxLayout, \
    QSpacerItem, QSizePolicy

from pyspinw.gui.helperwidgets.dockwidget import SpinWDockWidget
from pyspinw.gui.symmetry import SymmetryWidget
from pyspinw.gui.unitcell import UnitCellWidget


class LatticeParameters(SpinWDockWidget):
    def __init__(self, parent=None):

        super().__init__(parent=parent)

        self.setWindowTitle("Lattice")

        main_widget = QWidget()

        layout = QVBoxLayout()
        main_widget.setLayout(layout)

        self.symmetry_widget = SymmetryWidget(parent=self)
        self.unit_cell_widget = UnitCellWidget(parent=self)

        # Slots and signals

        self.symmetry_widget.symmetry_changed.connect(self._on_crystal_system_changed)
        self.unit_cell_widget.auto_update_lattice_system_request.connect(self.symmetry_widget.lattice_autoset)

        # Layout

        layout.addWidget(QLabel("Symmetry"))
        layout.addWidget(self.symmetry_widget)

        layout.addWidget(QLabel("Unit Cell"))
        layout.addWidget(self.unit_cell_widget)

        layout.addSpacerItem(QSpacerItem(0, 0, QSizePolicy.Minimum, QSizePolicy.Expanding))

        self.setWidget(main_widget)

    def _on_crystal_system_changed(self):
        """ Called when the crystal system is changed, updates the visible boxes"""
        self.unit_cell_widget.crystal_system = self.symmetry_widget.current_lattice_system

if __name__ == "__main__":
    app = QApplication([])

    parameter_widget = LatticeParameters()
    parameter_widget.show()

    app.exec_()
