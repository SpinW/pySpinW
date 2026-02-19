""" Text based display of the system """

from collections import defaultdict

import numpy as np
from PySide6.QtCore import Signal, QItemSelection, QItemSelectionModel
from PySide6.QtGui import QStandardItemModel, Qt, QStandardItem, QColor, QFont, QGuiApplication
from PySide6.QtWidgets import QWidget, QSplitter, QTextEdit, QTreeView, QVBoxLayout

from pyspinw.anisotropy import Anisotropy
from pyspinw.gui.render_model import RenderModel, RenderCoupling


class DisplayItem(QStandardItem):
    """ Slightly modified display item"""
    def __init__(self, text: str, ids: list[int]):
        super().__init__(text)
        self.ids = ids


def format_triple(triple):
    """ Rendering for tripples"""
    a,b,c = triple
    return f"({a:.3g}, {b:.3g}, {c:.3g})"

def format_moment_data(moment_data):
    """ Rendering for moment data"""
    if moment_data.shape[0] == 1:
        return format_triple(moment_data[0, :])
    else:
        return "--"

def format_anisotropies(anisotropies: list[Anisotropy]):
    if anisotropies:
        return "; ".join([a.parameter_string for a in anisotropies])
    else:
        return "None"

def format_g(g_matrix: np.ndarray):
    """ Formatting for g-factors"""
    diag = np.diagonal(g_matrix)
    if np.all(g_matrix == np.diag(diag)):
        a,b,c = diag

        if a==b and b==c:
            return f"{a:.3g}"
        else:
            return f"{a:.3g}, {b:.3g}, {c:.3g}"
    else:
        m = g_matrix.reshape(-1)
        return (f"[[{m[0]:.3g}, {m[1]:.3g}, {m[2]:.3g}],"
                f" [{m[3]:.3g}, {m[4]:.3g}, {m[5]:.3g}],"
                f" [{m[6]:.3g}, {m[7]:.3g}, {m[8]:.3g}]]")

class ParameterTable(QTreeView):
    """ Derivative class for both the text based site viewer and text based coupling viewer """

    hoverChanged = Signal()

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.entered.connect(self.on_item_hovered)

        self.hover_ids = []

        self.setMouseTracking(True)
        self.setSelectionBehavior(QTreeView.SelectRows)
        self.setSelectionMode(QTreeView.SelectionMode.ExtendedSelection)

        self.setStyleSheet("""
            QTreeView::item:selected {
                background-color: #DD9911;
                color: black;
            } """)


    def on_item_hovered(self, index):
        """ Called when item is hovered over"""
        if index.isValid():

            self.hover_ids = self.model().itemFromIndex(index).ids

            self.hoverChanged.emit()

    def leaveEvent(self, event):
        """ Qt override, mouse leaves widget"""

        self.hover_ids = []
        self.hoverChanged.emit()

        super().leaveEvent(event)



class TextDisplay(QWidget):
    """ Display of items in text form """

    hoverChanged = Signal()
    selectionChanged = Signal()

    _bold_font = QFont()
    _bold_font.setBold(True)

    _normal_font = QFont()

    def __init__(self, render_model: RenderModel, parent=None):

        super().__init__(parent=parent)


        self.default_color = QColor()
        self.hover_color = QColor("#ffcc00")
        self.selected_color = QColor()
        self.selected_hover_color = QColor()


        #
        # Sites and anisotropies
        #

        # mapping from site UID to anisotopies
        site_uid_to_anisotropies = defaultdict(list)

        for anisotropy in render_model.expanded.anisotropies:
            site_uid_to_anisotropies[anisotropy.site.unique_id].append(anisotropy)

        for anisotropy in render_model.hamiltonian.anisotropies:
            site_uid_to_anisotropies[anisotropy.site.unique_id].append(anisotropy)

        # Main set up
        self.site_tree = ParameterTable()

        site_model = QStandardItemModel()
        site_model.setHorizontalHeaderLabels([
            "Name", "Lattice Pos.", "Cartesian Pos.", "Moment", "Anisotropies", "g-Factor"])

        site_root = site_model.invisibleRootItem()

        small_cell = render_model.hamiltonian.structure.unit_cell
        big_cell = render_model.expanded.structure.unit_cell

        self.site_render_id_to_index = defaultdict(list)
        self.site_index_to_row_item: dict[int, DisplayItem] = {}

        #
        # Fill out site data
        #

        index = 0 # We need our own flat indexing
        for original_index, site in enumerate(render_model.original_sites):

            parent_ids = []

            unexpanded_name = DisplayItem(site.name, parent_ids)
            unexpanded_pos = DisplayItem(format_triple(site.ijk), parent_ids)
            unexpanded_cart = DisplayItem(format_triple(small_cell.fractional_to_cartesian(site.ijk)), parent_ids)
            unexpanded_moment = DisplayItem(format_moment_data(site.moment_data), parent_ids)
            unexpanded_anisotropies = DisplayItem(
                format_anisotropies(site_uid_to_anisotropies[site.unique_id]), parent_ids)
            unexpanded_gfactor = DisplayItem(format_g(site.g), parent_ids)

            parent_index = index

            self.site_index_to_row_item[index] = unexpanded_name

            index += 1

            for expanded_index in render_model.original_index_to_expanded[original_index]:

                expanded_render_site = render_model.sites[expanded_index]
                expanded_site = expanded_render_site.site

                ids=[render_model.sites[expanded_index].render_id]
                parent_ids += ids # This doesn't replace the list, adds to existing one

                expanded_name = DisplayItem(f"{expanded_site.name} {expanded_render_site.offset}", ids=ids)
                expanded_pos = DisplayItem(format_triple(expanded_site.ijk), ids=ids)
                expanded_cart = DisplayItem(format_triple(big_cell.fractional_to_cartesian(expanded_site.ijk)), ids=ids)
                expanded_moment = DisplayItem(format_triple(expanded_site.base_moment), ids=ids)
                expanded_anisotropies = DisplayItem(
                    format_anisotropies(site_uid_to_anisotropies[expanded_site.unique_id]), ids=ids)
                expanded_gfactor = DisplayItem(format_g(expanded_site.g), ids=ids)

                unexpanded_name.appendRow([
                    expanded_name,
                    expanded_pos,
                    expanded_cart,
                    expanded_moment,
                    expanded_anisotropies,
                    expanded_gfactor])

                self.site_index_to_row_item[index] = expanded_name

                for id in ids:
                    self.site_render_id_to_index[id].append(parent_index)
                    self.site_render_id_to_index[id].append(index)

                index += 1


            site_root.appendRow([
                unexpanded_name,
                unexpanded_pos,
                unexpanded_cart,
                unexpanded_moment,
                unexpanded_anisotropies,
                unexpanded_gfactor])

        self.site_tree.setModel(site_model)


        #
        # Couplings
        #

        # Get mapping from original to expanded
        self.couplings_original_to_expanded = defaultdict(list)
        for expanded_index, original_index in enumerate(render_model.coupling_mapping):
            self.couplings_original_to_expanded[original_index].append(expanded_index)


        # Fill out the coupling table
        self.coupling_tree = ParameterTable()

        coupling_model = QStandardItemModel()
        coupling_model.setHorizontalHeaderLabels(["Name", "Offset", "Type", "Site 1", "Site 2", "Parameters"])

        coupling_root = coupling_model.invisibleRootItem()

        for original_index, original_coupling in enumerate(render_model.hamiltonian.couplings):
            parent_ids = []
            unexpanded_name = DisplayItem(original_coupling.name, parent_ids)

            for expanded_index in self.couplings_original_to_expanded[original_index]:
                expanded_coupling: RenderCoupling = render_model.couplings[expanded_index]

                # Set the render ids
                ids = [expanded_coupling.render_id]
                parent_ids += ids # This doesn't replace the list, adds to existing one

                expanded_name = DisplayItem(expanded_coupling.coupling.name, ids)

                unexpanded_name.appendRow([
                    expanded_name
                ])


            coupling_root.appendRow([
                unexpanded_name
            ])

        self.coupling_tree.setModel(coupling_model)
        self.coupling_tree.expandAll()

        #
        # Layout components
        #

        splitter = QSplitter(Qt.Vertical)
        splitter.addWidget(self.site_tree)
        splitter.addWidget(self.coupling_tree)

        layout = QVBoxLayout()
        layout.addWidget(splitter)
        layout.setSpacing(0)
        layout.setContentsMargins(0,0,0,0)

        self.setLayout(layout)

        #
        # Hovering stuff
        #

        self.hover_ids = []
        self.current_selection = []
        self.site_tree.hoverChanged.connect(self.on_hover_changed)
        self.site_tree.selectionModel().selectionChanged.connect(self.on_selection_changed)
        self.coupling_tree.selectionModel().selectionChanged.connect(self.on_selection_changed)

    def on_hover_changed(self):
        self.hover_ids = self.site_tree.hover_ids # + other one
        self.hoverChanged.emit()

    def on_selection_changed(self):

        render_ids = []

        model = self.site_tree.model()
        for index in self.site_tree.selectionModel().selectedIndexes():
            if index.column() == 0:
                render_ids += model.itemFromIndex(index).ids

        model = self.coupling_tree.model()
        for index in self.coupling_tree.selectionModel().selectedIndexes():
            if index.column() == 0:
                render_ids += model.itemFromIndex(index).ids


        render_ids = list(set(render_ids))

        self.current_selection = render_ids
        self.selectionChanged.emit()


    @staticmethod
    def _get_row_items(item: QStandardItem):
        parent = item.parent()  # None if top-level
        row = item.row()

        if parent is None:
            model = item.model()
            column_count = model.columnCount()
            return [model.item(row, col) for col in range(column_count)]
        else:
            column_count = parent.columnCount()
            return [parent.child(row, col) for col in range(column_count)]

    def set_row_boldness(self, item: QStandardItem, is_bold: bool):

        for item in TextDisplay._get_row_items(item):

            if is_bold:
                item.setData(self._bold_font, Qt.FontRole)
            else:
                item.setData(self._normal_font, Qt.FontRole)


    def set_hover(self, render_ids):

        for item in self.site_index_to_row_item.values():
            self.set_row_boldness(item, False)

        for render_id in render_ids:
            for index in self.site_render_id_to_index[render_id]:
                # is_child = self.site_index_to_is_child[index]
                item = self.site_index_to_row_item[index]
                self.set_row_boldness(item, True)

    def set_selection(self, render_ids: list[int]):

        # print("Set selection:", render_ids)

        # Get the items to select
        site_items = []
        for render_id in render_ids:
            indices = self.site_render_id_to_index[render_id] # Defaultdict, so no item will be empty list
            for index in indices:
                item = self.site_index_to_row_item[index]
                parent = item.parent()
                if parent is not None:
                    site_items.append(item)
                    self.site_tree.expand(parent.index())

        selection = QItemSelection()

        for item in site_items:
            index = item.index()
            selection.select(index, index)


        selection_model = self.site_tree.selectionModel()

        selection_model.blockSignals(True) # Block the signals when sending
        try:
            selection_model.select(selection, QItemSelectionModel.ClearAndSelect | QItemSelectionModel.Rows)
        finally:
            selection_model.blockSignals(False)

        # Repaint signals will have got lost, so we need to manually refresh
        self.site_tree.viewport().update()
