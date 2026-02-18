from collections import defaultdict

from PySide6.QtCore import Signal
from PySide6.QtGui import QStandardItemModel, Qt, QStandardItem, QColor, QFont
from PySide6.QtWidgets import QWidget, QSplitter, QTextEdit, QTreeView, QVBoxLayout

from pyspinw.gui.render_model import RenderModel

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

class ParameterTable(QTreeView):

    hoverChanged = Signal()

    def __init__(self, parent=None):
        super().__init__(parent=parent)
        self.entered.connect(self.on_item_hovered)

        self.hover_ids = []

        self.setMouseTracking(True)

    def on_item_hovered(self, index):
        if index.isValid():

            self.hover_ids = self.model().itemFromIndex(index).ids

            self.hoverChanged.emit()

    def leaveEvent(self, event):

        self.hover_ids = []
        self.hoverChanged.emit()

        super().leaveEvent(event)



class TextDisplay(QWidget):

    hoverChanged = Signal()

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

        self.site_tree = ParameterTable()

        site_model = QStandardItemModel()
        site_model.setHorizontalHeaderLabels(["Name", "Lattice Pos.", "Cartesian Pos.", "Moment"])

        site_root = site_model.invisibleRootItem()

        small_cell = render_model.hamiltonian.structure.unit_cell
        big_cell = render_model.expanded.structure.unit_cell

        self.site_render_id_to_index = defaultdict(list)
        self.site_index_to_row_item: dict[int, DisplayItem] = {}

        index = 0 # We need our own flat indexing
        for original_index, site in enumerate(render_model.original_sites):

            parent_ids = []

            unexpanded_name = DisplayItem(site.name, parent_ids)
            unexpanded_pos = DisplayItem(format_triple(site.ijk), parent_ids)
            unexpanded_cart = DisplayItem(format_triple(small_cell.fractional_to_cartesian(site.ijk)), parent_ids)
            unexpanded_moment = DisplayItem(format_moment_data(site.moment_data), parent_ids)

            parent_index = index

            self.site_index_to_row_item[index] = unexpanded_name

            index += 1

            for expanded_index in render_model.original_index_to_expanded[original_index]:

                expanded_site = render_model.sites[expanded_index].site

                ids=[render_model.sites[expanded_index].render_id]
                parent_ids += ids

                expanded_name = DisplayItem(expanded_site.name, ids=ids)
                expanded_pos = DisplayItem(format_triple(expanded_site.ijk), ids=ids)
                expanded_cart = DisplayItem(format_triple(big_cell.fractional_to_cartesian(expanded_site.ijk)), ids=ids)
                expanded_moment = DisplayItem(format_triple(expanded_site.base_moment), ids=ids)

                unexpanded_name.appendRow([expanded_name, expanded_pos, expanded_cart, expanded_moment])

                self.site_index_to_row_item[index] = expanded_name

                for id in ids:
                    self.site_render_id_to_index[id].append(parent_index)
                    self.site_render_id_to_index[id].append(index)

                index += 1


            site_root.appendRow([unexpanded_name, unexpanded_pos, unexpanded_cart, unexpanded_moment])

        print(self.site_render_id_to_index)
        print(self.site_index_to_row_item)

        self.site_tree.setSelectionMode(QTreeView.SelectionMode.ExtendedSelection)
        self.site_tree.setModel(site_model)
        # site_tree.expandAll()

        #
        # Couplings
        #

        coupling_tree = QTreeView()

        coupling_model = QStandardItemModel()
        coupling_model.setHorizontalHeaderLabels(["Name", "Type", "Site 1", "Site 2"])

        root = coupling_model.invisibleRootItem()

        animals = QStandardItem("Animals")

        mammals = QStandardItem("Mammals")
        mammals.appendRow(QStandardItem("Dog"))
        mammals.appendRow(QStandardItem("Cat"))

        birds = QStandardItem("Birds")
        birds.appendRow(QStandardItem("Eagle"))

        animals.appendRow(mammals)
        animals.appendRow(birds)

        root.appendRow(animals)

        coupling_tree.setModel(coupling_model)
        coupling_tree.expandAll()


        splitter = QSplitter(Qt.Vertical)


        splitter.addWidget(self.site_tree)
        splitter.addWidget(coupling_tree)

        layout = QVBoxLayout()
        layout.addWidget(splitter)
        layout.setSpacing(0)
        layout.setContentsMargins(0,0,0,0)

        self.setLayout(layout)

        #
        # Hovering stuff
        #

        self.hover_ids = []
        self.site_tree.hoverChanged.connect(self.on_hover_changed)

    def on_hover_changed(self):
        self.hover_ids = self.site_tree.hover_ids # + other one


        self.hoverChanged.emit()

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