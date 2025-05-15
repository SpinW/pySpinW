import numpy as np
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QWidget, QVBoxLayout, QGridLayout, QApplication, QLabel

from pyspinw.gui.helperwidgets.misc import QRightLabel, QLeftLabel
from pyspinw.gui.helperwidgets.numbers import FloatField, FloatField
from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell


class CoordinateEditor(QWidget):
    """ Widget for editing a site """

    changed = Signal()

    def __init__(self, unit_cell: UnitCell | None = None, site: LatticeSite | None = None, parent=None):
        super().__init__(parent=parent)

        self._unit_cell = UnitCell(1,1,1) if unit_cell is None else unit_cell
        self._site = site

        self._position_editable = False
        self._moment_editable = False

        # Fractional part

        self.fractional_widget = QWidget(parent=self)

        fractional_layout = QGridLayout()

        self.i_editor = FloatField(0, slider_bottom=0, slider_top=1)
        self.j_editor = FloatField(0, slider_bottom=0, slider_top=1)
        self.k_editor = FloatField(0, slider_bottom=0, slider_top=1)

        self.mi_editor = FloatField(0, slider_bottom=-10, slider_top=10)
        self.mj_editor = FloatField(0, slider_bottom=-10, slider_top=10)
        self.mk_editor = FloatField(0, slider_bottom=-10, slider_top=10)

        self.i_editor.changed.connect(self._on_fractional_or_cell_changed)
        self.j_editor.changed.connect(self._on_fractional_or_cell_changed)
        self.k_editor.changed.connect(self._on_fractional_or_cell_changed)

        self.mi_editor.changed.connect(self._on_fractional_or_cell_changed)
        self.mj_editor.changed.connect(self._on_fractional_or_cell_changed)
        self.mk_editor.changed.connect(self._on_fractional_or_cell_changed)

        fractional_layout.addWidget(QRightLabel("i"), 0, 0)
        fractional_layout.addWidget(QRightLabel("j"), 1, 0)
        fractional_layout.addWidget(QRightLabel("k"), 2, 0)

        fractional_layout.addWidget(self.i_editor, 0, 1)
        fractional_layout.addWidget(self.j_editor, 1, 1)
        fractional_layout.addWidget(self.k_editor, 2, 1)

        ## No units for ijk

        fractional_layout.addWidget(QWidget(minimumWidth=10), 0, 3)

        fractional_layout.addWidget(QRightLabel("m<sub>i</sub>"), 0, 4)
        fractional_layout.addWidget(QRightLabel("m<sub>j</sub>"), 1, 4)
        fractional_layout.addWidget(QRightLabel("m<sub>k</sub>"), 2, 4)

        fractional_layout.addWidget(self.mi_editor, 0, 5)
        fractional_layout.addWidget(self.mj_editor, 1, 5)
        fractional_layout.addWidget(self.mk_editor, 2, 5)

        fractional_layout.addWidget(QLeftLabel("μ<sub>B</sub>/Å"), 0, 6)
        fractional_layout.addWidget(QLeftLabel("μ<sub>B</sub>/Å"), 1, 6)
        fractional_layout.addWidget(QLeftLabel("μ<sub>B</sub>/Å"), 2, 6)

        self.fractional_widget.setLayout(fractional_layout)

        # Cartesian part

        self.x_editor = FloatField(0, allow_slider=False)
        self.y_editor = FloatField(0, allow_slider=False)
        self.z_editor = FloatField(0, allow_slider=False)

        self.mx_editor = FloatField(0, allow_slider=False)
        self.my_editor = FloatField(0, allow_slider=False)
        self.mz_editor = FloatField(0, allow_slider=False)

        self.x_editor.changed.connect(self._on_cartesian_changed)
        self.y_editor.changed.connect(self._on_cartesian_changed)
        self.z_editor.changed.connect(self._on_cartesian_changed)

        self.mx_editor.changed.connect(self._on_cartesian_changed)
        self.my_editor.changed.connect(self._on_cartesian_changed)
        self.mz_editor.changed.connect(self._on_cartesian_changed)

        self.cartesian_widget = QWidget(parent=self)

        cartesian_layout = QGridLayout()


        cartesian_layout.addWidget(QRightLabel("x"), 0, 0)
        cartesian_layout.addWidget(QRightLabel("y"), 1, 0)
        cartesian_layout.addWidget(QRightLabel("z"), 2, 0)

        cartesian_layout.addWidget(self.x_editor, 0, 1)
        cartesian_layout.addWidget(self.y_editor, 1, 1)
        cartesian_layout.addWidget(self.z_editor, 2, 1)


        cartesian_layout.addWidget(QLeftLabel("Å"), 0, 2)
        cartesian_layout.addWidget(QLeftLabel("Å"), 1, 2)
        cartesian_layout.addWidget(QLeftLabel("Å"), 2, 2)

        cartesian_layout.addWidget(QWidget(minimumWidth=10), 0, 3)

        cartesian_layout.addWidget(QRightLabel("m<sub>x</sub>"), 0, 4)
        cartesian_layout.addWidget(QRightLabel("m<sub>y</sub>"), 1, 4)
        cartesian_layout.addWidget(QRightLabel("m<sub>z</sub>"), 2, 4)

        cartesian_layout.addWidget(self.mx_editor, 0, 5)
        cartesian_layout.addWidget(self.my_editor, 1, 5)
        cartesian_layout.addWidget(self.mz_editor, 2, 5)

        cartesian_layout.addWidget(QLeftLabel("μ<sub>B</sub>"), 0, 6)
        cartesian_layout.addWidget(QLeftLabel("μ<sub>B</sub>"), 1, 6)
        cartesian_layout.addWidget(QLeftLabel("μ<sub>B</sub>"), 2, 6)


        self.cartesian_widget.setLayout(cartesian_layout)

        # Both parts

        main_layout = QVBoxLayout()
        main_layout.addWidget(QLabel("Fractional Coordinates"))
        main_layout.addWidget(self.fractional_widget)
        main_layout.addWidget(QLabel("Cartesian Coordinates"))
        main_layout.addWidget(self.cartesian_widget)

        self.setLayout(main_layout)

        # Update the xyz coordinates
        self.position_editable = False
        self.moment_editable = False
        self._on_fractional_or_cell_changed()


    @property
    def unit_cell(self) -> UnitCell:
        return self._unit_cell

    @unit_cell.setter
    def unit_cell(self, unit_cell):
        self._unit_cell = unit_cell
        self._on_fractional_or_cell_changed()

    @property
    def position_editable(self):
        return self._position_editable

    @position_editable.setter
    def position_editable(self, editable: bool):
        self._position_editable = editable

        self.i_editor.setEnabled(editable)
        self.j_editor.setEnabled(editable)
        self.k_editor.setEnabled(editable)

        self.x_editor.setEnabled(editable)
        self.y_editor.setEnabled(editable)
        self.z_editor.setEnabled(editable)

    @property
    def moment_editable(self):
        return self._moment_editable

    @moment_editable.setter
    def moment_editable(self, editable):
        self._moment_editable = editable

        self.mi_editor.setEnabled(editable)
        self.mj_editor.setEnabled(editable)
        self.mk_editor.setEnabled(editable)

        self.mx_editor.setEnabled(editable)
        self.my_editor.setEnabled(editable)
        self.mz_editor.setEnabled(editable)

    @property
    def fractional_position(self):
        return np.array([self.i_editor.value, self.j_editor.value, self.k_editor.value])

    @property
    def fractional_moment(self):
        return np.array([self.mi_editor.value, self.mj_editor.value, self.mk_editor.value])


    @property
    def cartesian_position(self):
        return np.array([self.x_editor.value, self.y_editor.value, self.z_editor.value])

    @property
    def cartesian_moment(self):
        return np.array([self.mx_editor.value, self.my_editor.value, self.mz_editor.value])

    @property
    def site(self):
        return LatticeSite(
            i=self.i_editor.value,
            j=self.j_editor.value,
            k=self.k_editor.value,
            mi=self.mi_editor.value,
            mj=self.mj_editor.value,
            mk=self.mk_editor.value)

    @site.setter
    def site(self, site: LatticeSite):
        self.i_editor.value = site.i
        self.j_editor.value = site.j
        self.k_editor.value = site.k

        self.mi_editor.value = site.mi
        self.mj_editor.value = site.mj
        self.mk_editor.value = site.mk

    def _on_fractional_or_cell_changed(self):

        # self.cartesian_widget.blockSignals(True)

        new_cartesian_position = self.unit_cell.fractional_to_cartesian(self.fractional_position)
        new_cartesian_moment = self.unit_cell.fractional_to_cartesian(self.fractional_moment)

        self.x_editor.value = new_cartesian_position[0]
        self.y_editor.value = new_cartesian_position[1]
        self.z_editor.value = new_cartesian_position[2]

        self.mx_editor.value = new_cartesian_moment[0]
        self.my_editor.value = new_cartesian_moment[1]
        self.mz_editor.value = new_cartesian_moment[2]

        # self.cartesian_widget.blockSignals(False)

        self.changed.emit()

    def _on_cartesian_changed(self):


        # self.fractional_widget.blockSignals(True)

        new_fractional_position = self.unit_cell.cartesian_to_fractional(self.cartesian_position)
        new_fractional_moment = self.unit_cell.cartesian_to_fractional(self.cartesian_moment)

        self.i_editor.value = new_fractional_position[0]
        self.j_editor.value = new_fractional_position[1]
        self.k_editor.value = new_fractional_position[2]

        self.mi_editor.value = new_fractional_moment[0]
        self.mj_editor.value = new_fractional_moment[1]
        self.mk_editor.value = new_fractional_moment[2]

        # self.fractional_widget.blockSignals(False)

        self.changed.emit()


if __name__ == "__main__":
    app = QApplication([])

    editor = CoordinateEditor()
    editor.unit_cell = UnitCell(2,2,2)
    editor.show()

    app.exec_()