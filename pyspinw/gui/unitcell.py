""" GUI components for unit cells"""

from PySide6.QtCore import Signal, Qt
from PySide6.QtWidgets import QWidget, QGridLayout, QLabel, QVBoxLayout

from pyspinw.gui.helperwidgets.floatfield import FloatField
from pyspinw.symmetry.system import lattice_systems, LatticeSystem, find_unit_cell_type
from pyspinw.symmetry.unitcell import UnitCell



class UnitCellWidget(QWidget):
    """ Unit cell editing widget"""

    unit_cell_changed = Signal(UnitCell)
    auto_update_lattice_system_request = Signal(LatticeSystem)

    def __init__(self, unit_cell: UnitCell | None = None, parent=None):
        super().__init__(parent)

        self._current_unit_cell = UnitCell(1,1,1, 90, 90, 90) if unit_cell is None else unit_cell
        self._crystal_system = lattice_systems[0]

        self.a = FloatField(self.current_unit_cell.a, bottom=0.0001, slider_top=50, parent=self)
        self.b = FloatField(self.current_unit_cell.b, bottom=0.0001, slider_top=50, parent=self)
        self.c = FloatField(self.current_unit_cell.c, bottom=0.0001, slider_top=50, parent=self)

        self.alpha = FloatField(self.current_unit_cell.alpha, 0, 180, parent=self)
        self.beta = FloatField(self.current_unit_cell.beta, 0, 180, parent=self)
        self.gamma = FloatField(self.current_unit_cell.gamma, 0, 180, parent=self)

        self.a.changed.connect(self._on_cell_or_system_changed)
        self.b.changed.connect(self._on_cell_or_system_changed)
        self.c.changed.connect(self._on_cell_or_system_changed)

        self.alpha.changed.connect(self._on_cell_or_system_changed)
        self.beta.changed.connect(self._on_cell_or_system_changed)
        self.gamma.changed.connect(self._on_cell_or_system_changed)


        outer_layout = QVBoxLayout()
        self.setLayout(outer_layout)

        values_widget = QWidget()
        outer_layout.addWidget(values_widget)

        layout = QGridLayout(parent=self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        values_widget.setLayout(layout)

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

        self.message_label = QLabel("")
        outer_layout.addWidget(self.message_label)
        self.message_label.setStyleSheet("color: red;")

        self.message_label.setWordWrap(True)

        self.message_label.setOpenExternalLinks(False)  # Don't open in browser
        self.message_label.setTextInteractionFlags(self.message_label.textInteractionFlags() |
                                           Qt.LinksAccessibleByMouse)

        self.message_label.linkActivated.connect(self._auto_update_lattice_system)

        self._on_cell_or_system_changed()

    @property
    def lattice_system(self) -> LatticeSystem:
        """ Get the current lattice system"""
        return self._crystal_system

    @lattice_system.setter
    def lattice_system(self, lattice_system: LatticeSystem):
        self._crystal_system = lattice_system

        free_parameters = lattice_system.free_parameters

        self._set_grid_row_visible(0, free_parameters.a)
        self._set_grid_row_visible(1, free_parameters.b)
        self._set_grid_row_visible(2, free_parameters.c)
        self._set_grid_row_visible(3, free_parameters.alpha)
        self._set_grid_row_visible(4, free_parameters.beta)
        self._set_grid_row_visible(5, free_parameters.gamma)

        self._on_cell_or_system_changed()

    @property
    def current_unit_cell(self) -> UnitCell:
        """ Getter for the current unit cell"""
        return self._current_unit_cell

    def _set_grid_row_visible(self, row: int, value: bool):
        for i in range(self.grid_layout.columnCount()):
            item = self.grid_layout.itemAtPosition(row, i)
            if item and item.widget():
                item.widget().setVisible(value)

    def _on_cell_or_system_changed(self):
        """ Called when any element of the unit cell is changed"""
        self._current_unit_cell = self.lattice_system.constrain(
                                    UnitCell(
                                        self.a.value,
                                        self.b.value,
                                        self.c.value,
                                        self.alpha.value,
                                        self.beta.value,
                                        self.gamma.value))

        self.message_label.setText(UnitCellWidget._validity_message(self.lattice_system, self._current_unit_cell))

        self.unit_cell_changed.emit(self.current_unit_cell)



    @staticmethod
    def _validity_message(suggested_type: LatticeSystem, unit_cell: UnitCell) -> str:

        violated_constraints = suggested_type.violated_negative_constraints(unit_cell)

        if len(violated_constraints) == 0:
            return ""

        else:
            left = violated_constraints[:-1]
            right = violated_constraints[-1]

            if len(left) > 0:
                constraint_string = ", ".join(left) + " and " + right

            else:
                constraint_string = right

            actual_cell = find_unit_cell_type(unit_cell)

            if len(actual_cell) == 0:
                message = f"Invalid unit cell parameters, expected {constraint_string}"

            else:

                message = (f"Cell parameters have higher symmetry than {suggested_type.name}, "
                           f"expected {constraint_string}.")

                if len(actual_cell) == 1:
                    message += f'It looks like {actual_cell[0].name}. '
                    message += f'<a href="{actual_cell[0].name}"><font color="orange">Click here to set.</font></a>'

            return message

    def _auto_update_lattice_system(self, link):
        """ Called when the auto update link is clicked """
        self.auto_update_lattice_system_request.emit(link)
