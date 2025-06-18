import sys
from dataclasses import dataclass

import numpy as np

from PySide6.QtCore import Qt, QModelIndex, QSize, Signal, QEvent, QObject
from PySide6.QtGui import QTextDocument, QFontMetrics, QDoubleValidator, QIcon
from PySide6.QtWidgets import QTableWidget, QApplication, QHeaderView, QStyleOptionHeader, QStyle, QStyleOptionViewItem, \
    QTableWidgetItem, QAbstractItemView, QStyledItemDelegate, QLineEdit, QWidget, QCheckBox, QHBoxLayout

from pyspinw.gui.decorated import DecoratedSite, InteractionFlags
from pyspinw.gui.helperwidgets.misc import QRightLabel
from pyspinw.gui.symmetry_settings import SymmetrySettings, DEFAULT_SYMMETRY
from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.tolerances import tolerances
from pyspinw.util import problematic_sites


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
    """ Create a table entry for numbers"""
    item = QTableWidgetItem()
    item.setData(Qt.ItemDataRole.DisplayRole, float(x))
    if not editable:
        item.setFlags(item.flags() & ~Qt.ItemIsEditable)
    return item

_implied_icon = QIcon.fromTheme("folder")  # or use QIcon("path/to/icon.png")

def implied_entry(implied: bool):
    """ Create a table entry that signifies whether the site is implied or not"""

    item = QTableWidgetItem()

    if implied:
        item.setIcon(_implied_icon)

    item.setFlags(item.flags() & ~Qt.ItemIsEditable)

    return item

def name_entry(name: str, editable: bool = True):
    """ Create a table entry for names of things """

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



class HoverEventFilter(QObject):
    """ Event filter for detecting hover events, signals the row hovered, or -1 if nothing hovered """

    def __init__(self, table: "SiteTable"):
        super().__init__()
        self.table = table

    def eventFilter(self, source, event):
        """ Event filter that picks out the mouse move over the table and interprets it"""
        if event.type() == QEvent.MouseMove and source is self.table.viewport():
            index = self.table.indexAt(event.pos())
            if index.isValid():
                self.table._on_hover(index.row()) # Not the nicest way, but signals don't work from here
            else:
                self.table._on_hover(-1)

        return super().eventFilter(source, event)

@dataclass
class SiteTableOptions:
    """ Table display options"""
    nonmagnetic: bool = True
    independent: bool = True
    implied: bool = True
    supercell: bool = True

class SiteTableViewOptions(QWidget):
    """ Widget to choose what is displayed on the table"""

    options_changed = Signal()

    def __init__(self, parent, options=SiteTableOptions()):
        super().__init__(parent)

        layout = QHBoxLayout()

        layout.addWidget(QRightLabel("Show:"))

        self.nonmagnetic = QCheckBox("Non-magnetic")
        self.nonmagnetic.setChecked(options.nonmagnetic)
        self.nonmagnetic.checkStateChanged.connect(self._on_checked)
        layout.addWidget(self.nonmagnetic)

        self.independent = QCheckBox("Independent")
        self.independent.setChecked(options.independent)
        self.independent.checkStateChanged.connect(self._on_checked)
        layout.addWidget(self.independent)

        self.implied = QCheckBox("Implied")
        self.implied.setChecked(options.implied)
        self.implied.checkStateChanged.connect(self._on_checked)
        layout.addWidget(self.implied)

        self.supercell = QCheckBox("Supercell")
        self.supercell.setChecked(options.supercell)
        self.supercell.checkStateChanged.connect(self._on_checked)
        layout.addWidget(self.supercell)


    @property
    def options(self):
        """ Currently selected view options"""
        return SiteTableOptions(
            nonmagnetic=self.nonmagnetic.isChecked(),
            independent=self.independent.isChecked(),
            implied=self.implied.isChecked(),
            supercell=self.supercell.isChecked())

    def _on_checked(self):
        """Called when check boxes are checked"""
        self.options_changed.emit()



class SiteTable(QTableWidget):
    """ Table that shows the lattice sites """

    graphics_relevant_change = Signal()
    magnetic_symmetry_broken = Signal() # Sent when the magnetic symmetry is broken
    magnetic_symmetry_ok = Signal() # Sent when there has been an update, and the magnetic symmetry is ok

    implied_selected_state_changed = Signal(bool) # Signals whether selection contains a mirror site
    non_implied_selected_state_changed = Signal(bool) # Signals whether selection is non-empty

    def __init__(self, symmetry: SymmetrySettings=DEFAULT_SYMMETRY, parent=None):

        super().__init__(parent=parent)

        self._symmetry = symmetry

        # Data for the sites, the implied sites, and references both ways for linking them
        self._sites: list[LatticeSite] = []
        self._implied_sites: list[LatticeSite] = []
        self._implied_site_to_site: list[int] = []
        self._site_to_implied_site: list[list[int]] = []

        self.setRowCount(0)
        self.setColumnCount(14)

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

        # Multi-Row selection behaviour
        self.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.setSelectionMode(QAbstractItemView.MultiSelection)

        # callback for item selection
        self.itemSelectionChanged.connect(self._on_select)


        # Deligates for dealing with editing, needs to be kept in a list or python will dispose them
        #  and bad, confusing things will happen
        self._column_deligates = [FloatValidatorDelegate() for i in range(14)]

        for i in range(2,14):
            self.setItemDelegateForColumn(i, self._column_deligates[i])

        self.verticalHeader().hide()

        self.itemChanged.connect(self._on_item_changed)

        # Hovering
        self.setMouseTracking(True)
        self._hover_filter = HoverEventFilter(self)
        self.viewport().installEventFilter(self._hover_filter)
        self._current_hover_row = -1

        # Scale the columns
        self.resizeColumnsToContents()

    @property
    def sites(self) -> list[LatticeSite]:
        return self._sites

    @sites.setter
    def sites(self, sites: list[LatticeSite]):
        self._sites = sites
        self._update_entries()

    def add_site(self, site: LatticeSite):
        """ Add a site to the table"""
        self._sites.append(site)
        self._update_entries()

    def remove_site(self, index):
        """ Remove a site at a specified index"""
        self._sites.pop(index)
        self._update_entries()

    @property
    def symmetry(self):
        """ Get the current SymmetrySettings object"""
        return self._symmetry

    @symmetry.setter
    def symmetry(self, symmetry: SymmetrySettings):
        """ Set the current SymmetrySettings object"""
        self._symmetry = symmetry
        self._update_entries()


    def _update_sites(self):
        """ Update the sites, this will recalculate the implied sites and the arrays that
        say how they are referenced to each other"""
        implied_sites = []
        implied_site_to_site = []
        site_to_implied_site = []

        count = 0
        for site_index, site in enumerate(self._sites):
            extra_sites = self.symmetry.magnetic_group.duplicates(site)

            end_count = count + len(extra_sites)

            implied_sites += extra_sites
            implied_site_to_site += [site_index for _ in extra_sites]

            site_to_implied_site.append([i for i in range(count, end_count)])

            count = end_count

        self._implied_sites = implied_sites
        self._implied_site_to_site = implied_site_to_site
        self._site_to_implied_site = site_to_implied_site


    def _magnetic_symmetry_broken(self) -> bool:
        """ Check to see if the sites are consistent with the magnetic symmetry

        If the symmetry is broken, we will find there are duplicate sites with opposite spins
        """

        return len(problematic_sites(self._sites, self._implied_sites, self._site_to_implied_site)) > 0

    def _magnetic_symmetry_check(self):
        """ This finds out if the symmetry is broken, and sends the appropriate signals"""

        if self._magnetic_symmetry_broken():
            self.magnetic_symmetry_broken.emit()
        else:
            self.magnetic_symmetry_ok.emit()

    def _update_entries(self):
        """ Update the table entries """
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

        # See if the magnetic symmetry is ok
        self._magnetic_symmetry_check()

        # Main place to call this for actual changes
        self.graphics_relevant_change.emit()

    def magnetic_symmetry_autofix(self):
        """ Fix problems with magnetic symmetry automatically """
        bad_site_indices = problematic_sites(self._sites, self._implied_sites, self._site_to_implied_site)
        new_sites = []
        for idx, site in enumerate(self._sites):
            if idx in bad_site_indices:
                new_sites.append(LatticeSite.create(site.i, site.j, site.k, 0, 0, 0, name=site.name))
            else:
                new_sites.append(site)

        # Setting the sites should update everything relevant
        self.sites = new_sites

    @property
    def unit_cell(self):
        return self._symmetry.unit_cell

    @unit_cell.setter
    def unit_cell(self, unit_cell):
        raise Exception("Can't set the unit cell like this, use .symmetry with a SymmetrySettings object")

    def _on_hover(self, row: int):
        if self._current_hover_row != row:
            self._current_hover_row = row
            self.graphics_relevant_change.emit()

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
            new_site = LatticeSite.create(
                ijk[0], ijk[1], ijk[2],
                mijk[0], mijk[1], mijk[2],
                name = name)

        elif 7 <= col < 10:

            ijk_from_xyz = self._unit_cell.cartesian_to_fractional(xyz)

            new_site = LatticeSite.create(
                ijk_from_xyz[0], ijk_from_xyz[1], ijk_from_xyz[2],
                mijk[0], mikj[1], mijk[2],
                name=name)

        elif 10 <= col < 13:

            mijk_from_mxyz = self._unit_cell.cartesian_to_fractional(mxyz)

            new_site = LatticeSite.create(
                ijk[0], ijk[1], ijk[2],
                mijk_from_mxyz[0], mijk_from_mxyz[1], mijk_from_mxyz[2],
                name=name)

        else:
            return

        # replace the site in the list, and update

        self._sites.pop(row)
        self._sites.insert(row, new_site)


        self._update_entries()

    def _on_select(self):
        """ Event that happens on selection"""


        selected_indexes = set(idx.row() for idx in self.selectedIndexes())

        non_implied_selected = any([index < len(self._sites) for index in selected_indexes])
        self.non_implied_selected_state_changed.emit(non_implied_selected)

        inequivalent_selected = any([index >= len(self._sites) for index in selected_indexes])
        self.implied_selected_state_changed.emit(inequivalent_selected)


        # Inform graphics dependents
        self.graphics_relevant_change.emit()


    @property
    def _implicitness_marked_all_sites(self) -> list[LatticeSite, bool]:
        return [(site, False) for site in self._sites] + \
                [(site, True) for site in self._implied_sites]

    @property
    def selected_sites(self) -> list[LatticeSite]:
        """ Current selections """

        selected_indexes = sorted(list(set(idx.row() for idx in self.selectedIndexes())))

        sites = []
        n_real = len(self._sites)
        for index in selected_indexes:
            if index < n_real:
                sites.append(self._sites[index])
            else:
                sites.append(self._implied_sites[index - n_real])

        return sites


    @property
    def sites_for_drawing(self) -> list[DecoratedSite]:
        """ List of sites that will get sent to the graphics, these are decorated with information
        about the current selection"""

        out = []
        selected_indexes = set(idx.row() for idx in self.selectedIndexes())

        # Find which sites are in the group associated with the current hover
        if self._current_hover_row == -1:
            hover_indices = []
        else:

            n_sites = len(self._sites)

            # Find the base index
            if self._current_hover_row < n_sites:
                hover_index = self._current_hover_row
            else:
                hover_index = self._implied_site_to_site[self._current_hover_row - n_sites]

            # Create the list of associated sites
            hover_indices = [hover_index] + [idx + + n_sites for idx in self._site_to_implied_site[hover_index]]


        # Create the list of sites
        for i, (site, implied) in enumerate(self._implicitness_marked_all_sites):

            hover = self._current_hover_row == i
            selected = i in selected_indexes
            implied_hover = i in hover_indices

            # TODO: Add extra copies for edges of unit cells here, this will need to account for the extended cell,
            #       however, it might be better to do this elsewhere

            flags = InteractionFlags(
                hover=hover,
                selected=selected,
                implied=implied,
                implied_hover=implied_hover)

            # Get the position in the unit cell

            position = self.unit_cell.fractional_to_cartesian(site.ijk)
            moment = self.unit_cell.fractional_to_cartesian(site.m)

            # Add to the list
            out.append(DecoratedSite(
                site=site,
                position=position,
                moment=moment,
                flags=flags
            ))

        # print("render sites")
        # for site in out:
        #     print(site)

        return out

if __name__ == "__main__":
    app = QApplication([])

    site_table = SiteTable()

    site_table.unit_cell = UnitCell(2,3,4)

    site_table.add_site(LatticeSite(1,1,1))
    site_table.add_site(LatticeSite(1,1,2))
    site_table.add_site(LatticeSite(1,2,1))

    site_table.show()


    app.exec_()