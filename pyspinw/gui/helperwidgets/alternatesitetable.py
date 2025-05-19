import sys

import numpy as np

from PySide6.QtCore import Qt, QModelIndex, QSize, Signal
from PySide6.QtGui import QTextDocument, QFontMetrics, QDoubleValidator, QIcon
from PySide6.QtWidgets import QTableWidget, QApplication, QHeaderView, QStyleOptionHeader, QStyle, QStyleOptionViewItem, \
    QTableWidgetItem, QAbstractItemView, QStyledItemDelegate, QLineEdit

from pyspinw.gui.symmetry_settings import SymmetrySettings, DEFAULT_SYMMETRY
from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell


class HtmlHeader(QHeaderView):
    """ Header that can render html """

    def paintSection(self, painter, rect, logicalIndex):
        painter.save()

        # Draw background and borders
        opt = QStyleOptionHeader()
        self.initStyleOption(opt)
        opt.rect = rect
        opt.section = logicalIndex
        self.style().drawControl(QStyle.CE_Header, opt, painter, self)

        # Draw HTML text
        doc = QTextDocument()
        doc.setHtml(self.model().headerData(logicalIndex, self.orientation(), Qt.DisplayRole))
        doc.setTextWidth(rect.width())

        painter.translate(rect.topLeft())
        ctx = doc.documentLayout().PaintContext()
        doc.documentLayout().draw(painter, ctx)

        painter.restore()

    def sectionSizeFromContents(self, logicalIndex):
        """ Needed for the header to be resized"""

        # Return width required by HTML text
        doc = QTextDocument()
        html = self.model().headerData(logicalIndex, self.orientation(), Qt.DisplayRole)
        doc.setHtml(html)
        size = doc.size()
        # Slight padding for spacing
        return QSize(int(size.width()) + 10, int(size.height()) + 6)


def numeric_entry(x: float, editable: bool=True):
    """ How are we formatting numbers in the table"""
    item = QTableWidgetItem()
    item.setData(Qt.ItemDataRole.DisplayRole, float(x))
    if not editable:
        item.setFlags(item.flags() & ~Qt.ItemIsEditable)
    return item

_implied_icon = QIcon.fromTheme("folder")  # or use QIcon("path/to/icon.png")

def implied_entry(implied: bool):

    item = QTableWidgetItem()

    if implied:
        item.setIcon(_implied_icon)

    item.setFlags(item.flags() & ~Qt.ItemIsEditable)

    return item

def name_entry(name: str, editable: bool = True):
    item = QTableWidgetItem(name)

    if not editable:
        item.setFlags(item.flags() & ~Qt.ItemIsEditable)

    return item

class FloatValidatorDelegate(QStyledItemDelegate):
    """ Provides QLineEdits with float validators """

    def createEditor(self, parent, option, index):

        if not index.isValid():
            return None

        editor = QLineEdit(parent)
        # validator = QDoubleValidator(bottom=-sys.float_info.max, top=sys.float_info.min, decimals=100, parent=editor)
        validator = QDoubleValidator(bottom=-1e-100, top=1e100, decimals=100, parent=editor)
        validator.setNotation(QDoubleValidator.StandardNotation)
        editor.setValidator(validator)
        return editor



class SiteTable(QTableWidget):

    site_selected = Signal()

    def __init__(self, symmetry: SymmetrySettings=DEFAULT_SYMMETRY, parent=None):

        super().__init__(parent=parent)

        self._symmetry =symmetry
        self._sites: list[LatticeSite] = []
        self._implied_sites: list[LatticeSite] = []

        self.setRowCount(0)
        self.setColumnCount(13)

        header = HtmlHeader(Qt.Horizontal)
        self.setHorizontalHeader(header)
        self.setHorizontalHeaderLabels(["", "Name",
                                        "i", "j", "k",
                                        "m<sub>i</sub> (μ<sub>B</sub>/Å)",
                                        "m<sub>j</sub> (μ<sub>B</sub>/Å)",
                                        "m<sub>k</sub> (μ<sub>B</sub>/Å)",
                                        "x (Å)", "y (Å)", "z (Å)",
                                        "m<sub>x</sub> (μ<sub>B</sub>)",
                                        "m<sub>y</sub> (μ<sub>B</sub>)",
                                        "m<sub>z</sub> (μ<sub>B</sub>)"])

        self.setSelectionMode(QAbstractItemView.SingleSelection)
        self.itemSelectionChanged.connect(self._on_selection)

        # Deligates for dealing with editing, needs to be kept in a list or python will dispose them
        #  and bad, confusing things will happen
        self._column_deligates = [FloatValidatorDelegate() for i in range(14)]

        for i in range(2,14):
            self.setItemDelegateForColumn(i, self._column_deligates[i])

        self.verticalHeader().hide()

        self.itemChanged.connect(self._on_item_changed)
        self.resizeColumnsToContents()

    def add_site(self, site: LatticeSite):
        self._sites.append(site)
        self._update_entries()

    def _add_implied_site(self, site: LatticeSite):
        self._implied_sites.append(site)
        self._update_entries()

    def remove_site(self, index):
        self._sites.pop(index)
        self._update_entries()

    @property
    def symmetry(self):
        return self._symmetry

    @symmetry.setter
    def symmetry(self, symmetry: SymmetrySettings):
        self._symmetry = symmetry
        self._update_entries()

    def _on_selection(self):
        self.site_selected.emit()

    def _update_sites(self):
        implied_sites = []
        for site in self._sites:
            implied_sites += self.symmetry.magnetic_group.duplicates(site)
        self._implied_sites = implied_sites

    def _update_entries(self):
        self.blockSignals(True)

        self._update_sites()

        self.setRowCount(len(self._sites) + len(self._implied_sites))


        for i, site in enumerate(self._sites):
            self.setItem(i, 0, implied_entry(False))

            self.setItem(i, 1, name_entry(str(site.name)))

            self.setItem(i, 2, numeric_entry(site.i))
            self.setItem(i, 3, numeric_entry(site.j))
            self.setItem(i, 4, numeric_entry(site.k))

            self.setItem(i, 5, numeric_entry(site.mi))
            self.setItem(i, 6, numeric_entry(site.mj))
            self.setItem(i, 7, numeric_entry(site.mk))

            xyz = self.unit_cell.fractional_to_cartesian(site.ijk)
            mxyz = self.unit_cell.fractional_to_cartesian(site.m)

            self.setItem(i, 8, numeric_entry(xyz[0]))
            self.setItem(i, 9, numeric_entry(xyz[1]))
            self.setItem(i, 10, numeric_entry(xyz[2]))

            self.setItem(i, 11, numeric_entry(mxyz[0]))
            self.setItem(i, 12, numeric_entry(mxyz[1]))
            self.setItem(i, 13, numeric_entry(mxyz[2]))

        for j, site in enumerate(self._implied_sites):
            i = len(self._sites) + j

            self.setItem(i, 0, implied_entry(True))

            self.setItem(i, 1, name_entry(str(site.name), editable=False))


            self.setItem(i, 2, numeric_entry(site.i, editable=False))
            self.setItem(i, 3, numeric_entry(site.j, editable=False))
            self.setItem(i, 4, numeric_entry(site.k, editable=False))

            self.setItem(i, 5, numeric_entry(site.mi, editable=False))
            self.setItem(i, 6, numeric_entry(site.mj, editable=False))
            self.setItem(i, 7, numeric_entry(site.mk, editable=False))

            xyz = self.unit_cell.fractional_to_cartesian(site.ijk)
            mxyz = self.unit_cell.fractional_to_cartesian(site.m)

            self.setItem(i, 8, numeric_entry(xyz[0], editable=False))
            self.setItem(i, 9, numeric_entry(xyz[1], editable=False))
            self.setItem(i, 10, numeric_entry(xyz[2], editable=False))

            self.setItem(i, 11, numeric_entry(mxyz[0], editable=False))
            self.setItem(i, 12, numeric_entry(mxyz[1], editable=False))
            self.setItem(i, 13, numeric_entry(mxyz[2], editable=False))


        self.resizeColumnsToContents()
        self.blockSignals(False)

    @property
    def unit_cell(self):
        return self._symmetry.unit_cell

    @unit_cell.setter
    def unit_cell(self, unit_cell):
        raise Exception("Can't set the unit cell like this, use .symmetry with a SymmetrySettings object")

    def _on_item_changed(self, item: QTableWidgetItem):

        row = item.row()
        col = item.column()

        name = self.item(row, 1).text()

        ijk = np.array([
            float(self.item(row, 2).text()),
            float(self.item(row, 3).text()),
            float(self.item(row, 4).text())
        ])

        mijk = np.array([
            float(self.item(row, 5).text()),
            float(self.item(row, 6).text()),
            float(self.item(row, 7).text())
        ])

        xyz = np.array([
            float(self.item(row, 8).text()),
            float(self.item(row, 9).text()),
            float(self.item(row, 10).text())
        ])

        mxyz = np.array([
            float(self.item(row, 11).text()),
            float(self.item(row, 12).text()),
            float(self.item(row, 13).text())
        ])

        new_site = None

        if col < 7:
            new_site = LatticeSite(
                ijk[0], ijk[1], ijk[2],
                mijk[0], mijk[1], mijk[2],
                name = name)

        elif 7 <= col < 10:

            ijk_from_xyz = self._unit_cell.cartesian_to_fractional(xyz)

            new_site = LatticeSite(
                ijk_from_xyz[0], ijk_from_xyz[1], ijk_from_xyz[2],
                mijk[0], mikj[1], mijk[2],
                name=name)

        elif 10 <= col < 13:

            mijk_from_mxyz = self._unit_cell.cartesian_to_fractional(mxyz)

            new_site = LatticeSite(
                ijk[0], ijk[1], ijk[2],
                mijk_from_mxyz[0], mijk_from_mxyz[1], mijk_from_mxyz[2],
                name=name)

        else:
            return

        # replace the site in the list, and update

        self._sites.pop(row)
        self._sites.insert(row, new_site)


        self._update_entries()

if __name__ == "__main__":
    app = QApplication([])

    site_table = SiteTable()

    site_table.unit_cell = UnitCell(2,3,4)

    site_table.add_site(LatticeSite(1,1,1))
    site_table.add_site(LatticeSite(1,1,2))
    site_table.add_site(LatticeSite(1,2,1))
    site_table._add_implied_site(LatticeSite(1,2,1))

    site_table.show()


    app.exec_()