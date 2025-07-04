""" Table to display couplings """

import numpy as np
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

def _vector_format(vector: np.ndarray):
    """ Format a vector for display"""
    return ", ".join([f"{element:.4g}" for element in vector])


class CouplingTable(QTableWidget):
    """ Table of couplings """

    def __init__(self, parent=None, editable=True):
        super().__init__(parent=parent)

        self._couplings: list[Coupling] = []

        self.setRowCount(0)
        self.setColumnCount(7)

        self.verticalHeader().setVisible(False)
        self._set_headers()

        self._editable = editable


    def _set_headers(self):

        self.setHorizontalHeaderLabels(["Name",
                                        "Type",
                                        "Site 1",
                                        "Site 2",
                                        "Supercell",
                                        "Lattice Vector",
                                        "Parameters"])


    def update_entries(self):
        """ Update the table entries (view) """
        self.blockSignals(True)

        self.clear()

        self._set_headers()
        self.setRowCount(len(self._couplings))

        for i, coupling in enumerate(self._couplings):

            self.setItem(i, 0, string_entry(coupling.name, editable=self._editable))
            self.setItem(i, 1, string_entry(coupling.coupling_type, editable=False))
            self.setItem(i, 2, string_entry(coupling.site_1.name, editable=False))
            self.setItem(i, 3, string_entry(coupling.site_2.name, editable=False))
            self.setItem(i, 4, string_entry(str(coupling.cell_offset.as_tuple), editable=False))
            self.setItem(i, 5, string_entry(_vector_format(coupling.lattice_vector), editable=False))
            self.setItem(i, 6, string_entry(coupling.parameter_string, editable=False))


        self.resizeColumnsToContents()
        self.blockSignals(False)

    def add_coupling(self, coupling: Coupling):
        """ Add a coupling to the table"""
        self._couplings.append(coupling)

        self.update_entries()

    def add_couplings(self, couplings: list[Coupling]):
        """ Add a list of couplings to the table """
        self._couplings += couplings

        self.update_entries()

    @property
    def couplings(self):
        """ Get current couplings"""
        return self._couplings

    @couplings.setter
    def couplings(self, couplings: list[Coupling]):
        """ Set the current couplings"""
        self._couplings = couplings
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
