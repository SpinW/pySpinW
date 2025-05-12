from PySide6.QtCore import Signal
from PySide6.QtWidgets import QWidget, QGridLayout, QLabel

from pyspinw.gui.helperwidgets.numbers import FloatField
from pyspinw.symmetry.system import crystal_systems, CrystalSystem
from pyspinw.symmetry.unitcell import UnitCell

class UnitCellWidget(QWidget):

    unit_cell_changed = Signal(UnitCell)

    def __init__(self, unit_cell: UnitCell | None = None, parent=None):
        super().__init__(parent)

        self._current_unit_cell = UnitCell(1,1,1, 90, 90, 90) if unit_cell is None else unit_cell
        self._crystal_system = crystal_systems[0]

        self.a = FloatField(self.current_unit_cell.a, bottom=0.0001, slider_top=50, parent=self)
        self.b = FloatField(self.current_unit_cell.b, bottom=0.0001, slider_top=50, parent=self)
        self.c = FloatField(self.current_unit_cell.c, bottom=0.0001, slider_top=50, parent=self)

        self.alpha = FloatField(self.current_unit_cell.alpha, 0, 180, parent=self)
        self.beta = FloatField(self.current_unit_cell.beta, 0, 180, parent=self)
        self.gamma = FloatField(self.current_unit_cell.gamma, 0, 180, parent=self)

        self.a.changed.connect(self._on_cell_changed)
        self.b.changed.connect(self._on_cell_changed)
        self.c.changed.connect(self._on_cell_changed)

        self.alpha.changed.connect(self._on_cell_changed)
        self.beta.changed.connect(self._on_cell_changed)
        self.gamma.changed.connect(self._on_cell_changed)

        layout = QGridLayout(parent=self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        self.setLayout(layout)

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

        self.grid_layout = layout

    @property
    def crystal_system(self) -> CrystalSystem:
        return self._crystal_system

    @crystal_system.setter
    def crystal_system(self, crystal_system: CrystalSystem):
        self._crystal_system = crystal_system

        free_parameters = crystal_system.free_parameters

        self._set_grid_row_visible(0, free_parameters.a)
        self._set_grid_row_visible(1, free_parameters.b)
        self._set_grid_row_visible(2, free_parameters.c)
        self._set_grid_row_visible(3, free_parameters.alpha)
        self._set_grid_row_visible(4, free_parameters.beta)
        self._set_grid_row_visible(5, free_parameters.gamma)

        self._on_cell_changed()

    @property
    def current_unit_cell(self) -> UnitCell:
        return self._current_unit_cell

    def _set_grid_row_visible(self, row: int, value: bool):
        for i in range(self.grid_layout.columnCount()):
            item = self.grid_layout.itemAtPosition(row, i)
            if item and item.widget():
                item.widget().setVisible(value)

    def _on_cell_changed(self):
        """ Called when any element of the unit cell is changed"""
        self._current_unit_cell = self.crystal_system.constrain(
                                    UnitCell(
                                        self.a.value,
                                        self.b.value,
                                        self.c.value,
                                        self.alpha.value,
                                        self.beta.value,
                                        self.gamma.value))

        self.unit_cell_changed.emit(self.current_unit_cell)
