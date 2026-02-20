""" Text based display of the system """

from collections import defaultdict

import numpy as np
from PySide6.QtCore import Signal, QItemSelection, QItemSelectionModel
from PySide6.QtGui import QStandardItemModel, Qt, QStandardItem, QColor, QFont, QGuiApplication
from PySide6.QtWidgets import QWidget, QSplitter, QTextEdit, QTreeView, QVBoxLayout

from pyspinw.anisotropy import Anisotropy
from pyspinw.gui.rendermodel import RenderModel, RenderCoupling


class DisplayItem(QStandardItem):
    """ Slightly modified display item"""

    def __init__(self, text: str, ids: list[int]):
        super().__init__(text)
        self.ids = ids


def format_triple(triple):
    """ Rendering for numeric triples"""
    a,b,c = triple
    return f"({a:.3g}, {b:.3g}, {c:.3g})"

def format_moment_data(moment_data):
    """ Rendering for moment data"""
    if moment_data.shape[0] == 1:
        return format_triple(moment_data[0, :])
    else:
        return "--"

def format_anisotropies(anisotropies: list[Anisotropy]):
    """ Formatting for the anisotropies column"""
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

        # Mapping from render ID to first row item, use default dict to avoid key errors
        self.site_render_id_to_row_item: dict[int, list[DisplayItem]] = defaultdict(list)

        # Mapping from expanded site to offset
        site_expanded_uid_to_offset: dict[int, tuple[int, int, int]] = {}

        #
        # Fill out site data
        #

        index = 0 # We need our own flat indexing
        for unexpanded_index, site in enumerate(render_model.original_sites):

            # Note: unexpanded type is LatticeSite, expanded type is RenderSite

            parent_ids = []

            unexpanded_name = DisplayItem(site.name, parent_ids)
            unexpanded_pos = DisplayItem(format_triple(site.ijk), parent_ids)
            unexpanded_cart = DisplayItem(format_triple(small_cell.fractional_to_cartesian(site.ijk)), parent_ids)
            unexpanded_moment = DisplayItem(format_moment_data(site.moment_data), parent_ids)
            unexpanded_anisotropies = DisplayItem(
                format_anisotropies(site_uid_to_anisotropies[site.unique_id]), parent_ids)
            unexpanded_gfactor = DisplayItem(format_g(site.g), parent_ids)


            for expanded_index in render_model.original_index_to_expanded[unexpanded_index]:

                expanded_render_site = render_model.sites[expanded_index]
                expanded_site = expanded_render_site.site

                ids=[render_model.sites[expanded_index].render_id]
                parent_ids += ids # This doesn't replace the list, adds to existing one

                # Get the offset, save, then use for name
                offset = expanded_render_site.offset
                site_expanded_uid_to_offset[expanded_site.unique_id] = offset
                expanded_name = DisplayItem(f"{expanded_site.name} {offset}", ids=ids)

                # rest of items
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

                for id in ids:
                    self.site_render_id_to_row_item[id] = [unexpanded_name, expanded_name]



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
        # Mapping from render ID to first row item, use default dict to avoid key errors
        self.coupling_render_id_to_row_item: dict[int, list[DisplayItem]] = defaultdict(list)

        # Get mapping from original to expanded
        self.couplings_original_to_expanded = defaultdict(list)
        for expanded_index, unexpanded_index in enumerate(render_model.coupling_mapping):
            self.couplings_original_to_expanded[unexpanded_index].append(expanded_index)


        # Fill out the coupling table
        self.coupling_tree = ParameterTable()

        coupling_model = QStandardItemModel()
        coupling_model.setHorizontalHeaderLabels(["Name", "Offset", "Type", "Site 1", "Site 2", "Parameters"])

        coupling_root = coupling_model.invisibleRootItem()

        for unexpanded_index, unexpanded_coupling in enumerate(render_model.hamiltonian.couplings):

            # Note unexpanded is type Coupling, expanded is type RenderCoupling

            parent_ids = []
            unexpanded_name = DisplayItem(unexpanded_coupling.name, parent_ids)
            unexpanded_offset = DisplayItem(str(unexpanded_coupling.cell_offset.as_tuple), parent_ids)
            unexpanded_type = DisplayItem(unexpanded_coupling.coupling_type, parent_ids)
            unexpanded_site_1 = DisplayItem(unexpanded_coupling.site_1.name, parent_ids)
            unexpanded_site_2 = DisplayItem(unexpanded_coupling.site_2.name, parent_ids)
            unexpanded_parameters = DisplayItem(unexpanded_coupling.parameter_string, parent_ids)


            for expanded_index in self.couplings_original_to_expanded[unexpanded_index]:
                expanded_coupling: RenderCoupling = render_model.couplings[expanded_index]

                # Set the render ids
                ids = [expanded_coupling.render_id]
                parent_ids += ids # This doesn't replace the list, adds to existing one

                site_1_offset = site_expanded_uid_to_offset[expanded_coupling.coupling.site_1.unique_id]
                site_1_name = f"{expanded_coupling.coupling.site_1.name} {site_1_offset}"
                site_2_offset = site_expanded_uid_to_offset[expanded_coupling.coupling.site_2.unique_id]
                site_2_name = f"{expanded_coupling.coupling.site_2.name} {site_2_offset}"

                expanded_name = DisplayItem(expanded_coupling.coupling.name, ids)
                expanded_offset = DisplayItem(str(expanded_coupling.coupling.cell_offset.as_tuple), ids)
                expanded_type = DisplayItem(expanded_coupling.coupling.coupling_type, ids)
                expanded_site_1 = DisplayItem(site_1_name, ids)
                expanded_site_2 = DisplayItem(site_2_name, ids)
                expanded_parameters = DisplayItem(expanded_coupling.coupling.parameter_string, ids)

                unexpanded_name.appendRow([
                    expanded_name,
                    expanded_offset,
                    expanded_type,
                    expanded_site_1,
                    expanded_site_2,
                    expanded_parameters
                ])

                # Save column zero items
                for id in ids:
                    self.coupling_render_id_to_row_item[id] = [unexpanded_name, expanded_name]


            coupling_root.appendRow([
                unexpanded_name,
                unexpanded_offset,
                unexpanded_type,
                unexpanded_site_1,
                unexpanded_site_2,
                unexpanded_parameters
            ])

        self.coupling_tree.setModel(coupling_model)

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
        # Hovering and selection stuff
        #

        self.hover_ids = []
        self.current_selection = []
        self.site_tree.hoverChanged.connect(self.on_hover_changed)
        self.site_tree.selectionModel().selectionChanged.connect(self.on_selection_changed)
        self.coupling_tree.selectionModel().selectionChanged.connect(self.on_selection_changed)

    def on_hover_changed(self):
        """ Called when items are hovered """
        self.hover_ids = self.site_tree.hover_ids # + other one
        self.hoverChanged.emit()

    def on_selection_changed(self):
        """ Called when selection on either of the trees changes"""
        #
        # Build a list of render ids
        #

        render_ids = []

        # Sites
        model = self.site_tree.model()
        for index in self.site_tree.selectionModel().selectedIndexes():
            if index.column() == 0:
                render_ids += model.itemFromIndex(index).ids

        # Couplings
        model = self.coupling_tree.model()
        for index in self.coupling_tree.selectionModel().selectedIndexes():
            if index.column() == 0:
                render_ids += model.itemFromIndex(index).ids

        # Remove duplicates
        render_ids = list(set(render_ids))

        # Make accessible, then signal anything listening
        self.current_selection = render_ids
        self.selectionChanged.emit()


    @staticmethod
    def _get_row_items(item: QStandardItem):
        """ Get items in a row """
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
        """ Set the boldness of the font for a given row"""
        for item in TextDisplay._get_row_items(item):

            if is_bold:
                item.setData(self._bold_font, Qt.FontRole)
            else:
                item.setData(self._normal_font, Qt.FontRole)


    def set_hover(self, render_ids: list[int]):
        """ Set the hover state for specified render_ids"""
        for items in self.site_render_id_to_row_item.values():
            # TODO Be more efficient
            for item in items:
                self.set_row_boldness(item, False)

        for render_id in render_ids:
            for item in self.site_render_id_to_row_item[render_id]:
                self.set_row_boldness(item, True)

        for items in self.coupling_render_id_to_row_item.values():
            # TODO Be more efficient
            for item in items:
                self.set_row_boldness(item, False)

        for render_id in render_ids:
            for item in self.coupling_render_id_to_row_item[render_id]:
                self.set_row_boldness(item, True)

    def set_selection(self, render_ids: list[int]):
        """ Set the selection based on render IDs"""
        # print("Set selection:", render_ids)

        #
        # Sites
        #

        selected_items = [] # Items to select

        # Get the selection items, and expand where needed
        for render_id in render_ids:
            for item in self.site_render_id_to_row_item[render_id]:
                parent = item.parent()
                if parent is not None:
                    selected_items.append(item)
                    self.site_tree.expand(parent.index())

        # Make the selection object
        selection = QItemSelection()

        for item in selected_items:
            index = item.index()
            selection.select(index, index)

        # Do selection

        selection_model = self.site_tree.selectionModel()

        selection_model.blockSignals(True) # Block the signals when sending
        try:
            selection_model.select(selection, QItemSelectionModel.ClearAndSelect | QItemSelectionModel.Rows)
        finally:
            selection_model.blockSignals(False)

        # Repaint signals will have got lost, so we need to manually refresh
        self.site_tree.viewport().update()

        #
        # Couplings
        #

        selected_items = [] # items to be selected

        # Get the selection items, and expand where needed
        for render_id in render_ids:
            for item in self.coupling_render_id_to_row_item[render_id]:
                parent = item.parent()
                if parent is not None:
                    selected_items.append(item)
                    self.coupling_tree.expand(parent.index())

        selection = QItemSelection()

        for item in selected_items:
            index = item.index()
            selection.select(index, index)

        # Set selection on tree
        selection_model = self.coupling_tree.selectionModel()

        selection_model.blockSignals(True)  # Block the signals when sending
        try:
            selection_model.select(selection, QItemSelectionModel.ClearAndSelect | QItemSelectionModel.Rows)
        finally:
            selection_model.blockSignals(False)

        # Repaint signals will have got lost, so we need to manually refresh
        self.coupling_tree.viewport().update()
