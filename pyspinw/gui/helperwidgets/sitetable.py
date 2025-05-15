from PySide6.QtCore import Qt, QModelIndex, QSize, Signal
from PySide6.QtGui import QTextDocument, QFontMetrics
from PySide6.QtWidgets import QTableWidget, QApplication, QHeaderView, QStyleOptionHeader, QStyle, QStyleOptionViewItem, \
    QTableWidgetItem, QAbstractItemView

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

def numeric_entry(x: float):
    """ How are we formatting numbers in the table"""
    return QTableWidgetItem(f"{x:.4g}")


class SiteTable(QTableWidget):

    site_selected = Signal()

    def __init__(self, unit_cell: UnitCell | None = None, parent=None):

        super().__init__(parent=parent)

        self._unit_cell = UnitCell(1,1,1) if unit_cell is None else unit_cell
        self._sites: list[LatticeSite] = []

        self.setRowCount(0)
        self.setColumnCount(13)

        header = HtmlHeader(Qt.Horizontal)
        self.setHorizontalHeader(header)
        self.setHorizontalHeaderLabels(["Name", "i", "j", "k", "m<sub>i</sub>", "m<sub>j</sub>", "m<sub>k</sub>", "x", "y", "z", "m<sub>x</sub>", "m<sub>y</sub>", "m<sub>z</sub>"])

        self.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.setSelectionMode(QAbstractItemView.SingleSelection)
        self.itemSelectionChanged.connect(self._on_selection)

        self.resizeColumnsToContents()

    def add_site(self, site: LatticeSite):
        self._sites.append(site)
        self._update_entries()

    def remove_site(self, index):
        self._sites.pop(index)
        self._update_entries()

    def _on_selection(self):
        self.site_selected.emit()

    def _update_entries(self):
        self.setRowCount(len(self._sites))
        for i, site in enumerate(self._sites):
            self.setItem(i, 0, QTableWidgetItem(str(site.name)))

            self.setItem(i, 1, numeric_entry(site.i))
            self.setItem(i, 2, numeric_entry(site.j))
            self.setItem(i, 3, numeric_entry(site.k))

            self.setItem(i, 4, numeric_entry(site.mi))
            self.setItem(i, 5, numeric_entry(site.mj))
            self.setItem(i, 6, numeric_entry(site.mk))

            xyz = self._unit_cell.fractional_to_cartesian(site.ijk)
            mxyz = self._unit_cell.fractional_to_cartesian(site.m)

            self.setItem(i, 7, numeric_entry(xyz[0]))
            self.setItem(i, 8, numeric_entry(xyz[1]))
            self.setItem(i, 9, numeric_entry(xyz[2]))

            self.setItem(i, 10, numeric_entry(mxyz[0]))
            self.setItem(i, 11, numeric_entry(mxyz[1]))
            self.setItem(i, 12, numeric_entry(mxyz[2]))

        self.resizeColumnsToContents()


if __name__ == "__main__":
    app = QApplication([])

    site_table = SiteTable()

    site_table.add_site(LatticeSite(1,1,1))
    site_table.add_site(LatticeSite(1,1,2))
    site_table.add_site(LatticeSite(1,2,1))

    site_table.show()


    app.exec_()