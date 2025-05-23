from PySide6.QtCore import Qt
from PySide6.QtWidgets import QTableWidget, QTableWidgetItem, QApplication

from pyspinw.coupling import Coupling, DMCoupling
from pyspinw.gui.cell_offsets import CellOffset
from pyspinw.site import LatticeSite


def string_entry(value: str, editable: bool = True):
    """ Create a table entry for names of things """

    item = QTableWidgetItem(value)

    if not editable:
        item.setFlags(item.flags() & ~Qt.ItemIsEditable)

    return item

class CouplingTable(QTableWidget):
    def __init__(self, parent=None):
        super().__init__(parent=parent)

        self._couplings: list[Coupling] = []

        self.setRowCount(0)
        self.setColumnCount(5)

        self.setHorizontalHeaderLabels(["Name",
                                        "Site 1",
                                        "Site 2",
                                        "Supercell"
                                        "Type",
                                        "Parameters"])


    def update_entries(self):

        self.blockSignals(True)

        self.setRowCount(len(self._couplings))


        for i, coupling in enumerate(self._couplings):

            self.setItem(i, 0, string_entry(coupling.name))
            self.setItem(i, 1, string_entry(coupling.site_1.name, editable=False))
            self.setItem(i, 2, string_entry(coupling.site_2.name, editable=False))
            self.setItem(i, 3, string_entry(str(coupling.cell_offset), editable=False))
            self.setItem(i, 4, string_entry(coupling.coupling_type, editable=False))
            self.setItem(i, 5, string_entry(coupling.parameter_string, editable=False))


        self.resizeColumnsToContents()
        self.blockSignals(False)

    def add_coupling(self, coupling: Coupling):
        self._couplings.append(coupling)

        self.update_entries()

if __name__ == "__main__":
    app = QApplication([])

    coupling_table = CouplingTable()

    a = LatticeSite(1, 1, 1, name="a")
    b = LatticeSite(1, 1, 1, name="b")
    c = LatticeSite(1, 1, 1, name="c")

    couplings = [
        DMCoupling(name="DM1", site_1 = a, site_2 = b, d_x=1, d_y=1, d_z=1, cell_offset=CellOffset(i=0,j=0,k=0))
    ]

    for coupling in couplings:
        coupling_table.add_coupling(coupling)

    coupling_table.show()

    app.exec_()
