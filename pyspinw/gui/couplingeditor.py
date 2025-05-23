from PySide6.QtWidgets import QWidget, QGridLayout, QApplication, QLineEdit, QComboBox

from pyspinw.coupling import Coupling, DMCoupling
from pyspinw.gui.cell_offsets import CellOffset
from pyspinw.gui.helperwidgets.couplingtypecombo import CouplingTypeCombo
from pyspinw.gui.helperwidgets.floatfield import FloatField
from pyspinw.gui.helperwidgets.index_field import IndexField
from pyspinw.gui.helperwidgets.misc import QLeftLabel, QRightLabel
from pyspinw.site import LatticeSite


class CouplingEditor(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        self._coupling: Coupling | None = None

        self.name_field = QLineEdit()
        self.i_field = QLineEdit()

        self.main_layout = QGridLayout()

        self.setLayout(self.main_layout)

    def _clear(self):
        for i in range(self.main_layout.count()):
            self.main_layout.itemAt(i).widget().deleteLater()


    def _update_parameters(self):
        self._clear()

        if self._coupling is None:
            return

        self.main_layout.addWidget(QRightLabel("Name"), 0, 0)
        self.main_layout.addWidget(QLineEdit(self._coupling.name), 0, 1)

        self.main_layout.addWidget(QRightLabel("Site 1"), 1, 0)
        self.main_layout.addWidget(QComboBox(), 1, 1)

        self.main_layout.addWidget(QRightLabel("Site 2"), 2, 0)
        self.main_layout.addWidget(QComboBox(), 2, 1)

        self.main_layout.addWidget(QRightLabel("Supercell"), 3, 0)
        self.main_layout.addWidget(IndexField(self._coupling.cell_offset.as_tuple), 3, 1)

        self.main_layout.addWidget(QRightLabel("Type"))
        self.main_layout.addWidget(CouplingTypeCombo(self._coupling.coupling_type))

        for i, parameter in enumerate(self._coupling.parameters):
            field = FloatField(self._coupling.__dict__[parameter], slider_bottom=-50, slider_top=50, parent=self)
            self.main_layout.addWidget(QRightLabel(parameter), i+5, 0)
            self.main_layout.addWidget(field, i+5, 1)

    @property
    def coupling(self):
        return self._coupling

    @coupling.setter
    def coupling(self, coupling):
        self._coupling = coupling

        self._update_parameters()


if __name__ == "__main__":
    app = QApplication([])

    coupling_editor = CouplingEditor()

    a = LatticeSite(1, 1, 1, name="a")
    b = LatticeSite(1, 1, 1, name="b")

    coupling = DMCoupling(name="DM1", site_1 = a, site_2 = b, d_x=1, d_y=1, d_z=1, cell_offset=CellOffset(i=0,j=0,k=0))

    coupling_editor.coupling = coupling

    coupling_editor.show()

    app.exec_()
