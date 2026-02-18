from PySide6.QtCore import Signal
from PySide6.QtGui import QStandardItemModel, Qt, QStandardItem
from PySide6.QtWidgets import QWidget, QSplitter, QTextEdit, QTreeView, QVBoxLayout

from pyspinw.gui.render_model import RenderModel

class DisplayItem(QStandardItem):
    def __init__(self, text: str, ids: list[int]):
        super().__init__(text)
        self.ids = ids


def format_tripple(tripple):
    a,b,c = tripple
    return f"({a:.3g}, {b:.3g}, {c:.3g})"

def format_moment_data(moment_data):
    if moment_data.shape[0] == 1:
        return format_tripple(moment_data[0, :])
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

    def __init__(self, render_model: RenderModel, parent=None):

        super().__init__(parent=parent)

        #
        # Sites and anisotropies
        #

        self.site_tree = ParameterTable()

        site_model = QStandardItemModel()
        site_model.setHorizontalHeaderLabels(["Name", "Lattice Pos.", "Cartesian Pos.", "Moment"])

        site_root = site_model.invisibleRootItem()

        small_cell = render_model.hamiltonian.structure.unit_cell
        big_cell = render_model.expanded.structure.unit_cell

        for original_index, site in enumerate(render_model.original_sites):

            parent_ids = []

            unexpanded_name = DisplayItem(site.name, parent_ids)
            unexpanded_pos = DisplayItem(format_tripple(site.ijk), parent_ids)
            unexpanded_cart = DisplayItem(format_tripple(small_cell.fractional_to_cartesian(site.ijk)), parent_ids)
            unexpanded_moment = DisplayItem(format_moment_data(site.moment_data), parent_ids)

            for expanded_index in render_model.original_index_to_expanded[original_index]:

                expanded_site = render_model.sites[expanded_index].site

                ids=[render_model.sites[expanded_index].render_id]
                parent_ids += ids

                expanded_name = DisplayItem(expanded_site.name, ids=ids)
                expanded_pos = DisplayItem(format_tripple(expanded_site.ijk), ids=ids)
                expanded_cart = DisplayItem(format_tripple(big_cell.fractional_to_cartesian(expanded_site.ijk)), ids=ids)
                expanded_moment = DisplayItem(format_tripple(expanded_site.base_moment), ids=ids)

                unexpanded_name.appendRow([expanded_name, expanded_pos, expanded_cart, expanded_moment])

            site_root.appendRow([unexpanded_name, unexpanded_pos, unexpanded_cart, unexpanded_moment])

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