from PySide6.QtCore import Signal
from PySide6.QtWidgets import QBoxLayout, QGridLayout, QWidget, QApplication, QLabel, QDockWidget

from pyspinw.gui.helperwidgets.dockwidget import SpinWDockWidget
from pyspinw.gui.helperwidgets.numbers import FloatField
from pyspinw.unitcell import UnitCell


class LatticeParameters(SpinWDockWidget):

    unit_cell_changed = Signal(UnitCell)

    def __init__(self, unit_cell: UnitCell | None = None, parent=None):
        super().__init__(parent)

        self.setWindowTitle("Lattice")
        self.current_unit_cell = UnitCell(1,1,1, 90, 90, 90) if unit_cell is None else unit_cell

        self.a = FloatField(self.current_unit_cell.a, bottom=0.0001, slider_top=50)
        self.b = FloatField(self.current_unit_cell.b, bottom=0.0001, slider_top=50)
        self.c = FloatField(self.current_unit_cell.c, bottom=0.0001, slider_top=50)

        self.alpha = FloatField(self.current_unit_cell.alpha, 0, 180)
        self.beta = FloatField(self.current_unit_cell.beta, 0, 180)
        self.gamma = FloatField(self.current_unit_cell.gamma, 0, 180)

        self.a.changed.connect(self._on_lattice_changed)
        self.b.changed.connect(self._on_lattice_changed)
        self.c.changed.connect(self._on_lattice_changed)

        self.alpha.changed.connect(self._on_lattice_changed)
        self.beta.changed.connect(self._on_lattice_changed)
        self.gamma.changed.connect(self._on_lattice_changed)

        main_widget = QWidget()
        self.setWidget(main_widget)

        layout = QGridLayout(parent=self)
        main_widget.setLayout(layout)

        layout.addWidget(QLabel("a"), 0, 0)
        layout.addWidget(QLabel("b"), 1, 0)
        layout.addWidget(QLabel("c"), 2, 0)

        layout.addWidget(QLabel("α"), 3, 0)
        layout.addWidget(QLabel("β"), 4, 0)
        layout.addWidget(QLabel("γ"), 5, 0)

        layout.addWidget(self.a, 0, 1)
        layout.addWidget(self.b, 1, 1)
        layout.addWidget(self.c, 2, 1)

        layout.addWidget(self.alpha, 3, 1)
        layout.addWidget(self.beta, 4, 1)
        layout.addWidget(self.gamma, 5, 1)

        layout.addWidget(QLabel("Å"), 0, 2)
        layout.addWidget(QLabel("Å"), 1, 2)
        layout.addWidget(QLabel("Å"), 2, 2)

        layout.addWidget(QLabel("°"), 3, 2)
        layout.addWidget(QLabel("°"), 4, 2)
        layout.addWidget(QLabel("°"), 5, 2)

    def _on_lattice_changed(self):
        self.current_unit_cell = UnitCell(
                                    self.a.value,
                                    self.b.value,
                                    self.c.value,
                                    self.alpha.value,
                                    self.beta.value,
                                    self.gamma.value)

        self.unit_cell_changed.emit(self.current_unit_cell)

if __name__ == "__main__":
    app = QApplication([])

    parameter_widget = LatticeParameters()
    parameter_widget.show()

    app.exec_()
