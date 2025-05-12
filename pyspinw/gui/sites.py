from PySide6.QtWidgets import QWidget, QVBoxLayout, QTreeWidget, QLabel, QApplication, QTreeWidgetItem, \
    QAbstractItemView, QListWidget, QListWidgetItem

from pyspinw.gui.helperwidgets.dockwidget import SpinWDockWidget
from pyspinw.site import LatticeSite
from pyspinw.symmetry.group import SymmetryGroup


class AsymmetricSiteItem(QListWidgetItem):
    def __init__(self, site: LatticeSite):
        super().__init__("A site")

        self.site = site



class AsymmetricCellSites(QWidget):
    """ Sites that belong to the asymmetric cell """
    def __init__(self, parent=None):
        super().__init__(parent=parent)

        layout = QVBoxLayout()
        self.setLayout(layout)

        self.tree = QListWidget(parent=self)

        self.tree.setSelectionBehavior(QAbstractItemView.SelectItems)
        # self.tree.selectionModel().selectionChanged.connect(self.onSelected)
        # self.tree.setHeaderHidden(True)

        layout.addWidget(QLabel("Asymmetric Cell"))
        layout.addWidget(self.tree)


    def add_site(self, site: LatticeSite):
        self.tree.addItem(AsymmetricSiteItem(site))



class SymmetrySites(QWidget):
    """ Tree of sites that are implied by symmetry """
    def __init__(self, parent=None):
        super().__init__(parent=parent)

        layout = QVBoxLayout()
        self.setLayout(layout)

        self.tree = QTreeWidget(parent=self)
        self.tree.setHeaderHidden(True)

        layout.addWidget(QLabel("Implied Sites"))
        layout.addWidget(self.tree)

    def update(self, asymmetric_sites: list[LatticeSite], symmetry: SymmetryGroup):
        """ Update the implicit sites given by the symmetry """

        # Does this update affect what sites there are, rather than the locations, if so, we want to give a warning
        # We don't want to warn if we've deleted/added something from the initial list, just if its moved to a special point
        # Or if we have changed the symmetry
        # We should also warn if we're deleting something that affects a constraint


class SiteEditor(SpinWDockWidget):
    """ Editor dock window for magnetic sites"""
    def __init__(self, parent=None):
        super().__init__(parent=parent)

        self.setWindowTitle("Sites")

        widget = QWidget()
        layout = QVBoxLayout()

        self.setWidget(widget)
        widget.setLayout(layout)

        self.asymmetric_cell_sites = AsymmetricCellSites(self)
        self.symmetry_cell_sites = SymmetrySites(self)

        layout.addWidget(self.asymmetric_cell_sites)
        layout.addWidget(self.symmetry_cell_sites)


if __name__ == "__main__":
    app = QApplication([])

    site_editor = SiteEditor()
    site_editor.asymmetric_cell_sites.add_site(LatticeSite(0.5, 0.5, 0.5))
    site_editor.asymmetric_cell_sites.add_site(LatticeSite(0.5, 0.5, 0.5))
    site_editor.asymmetric_cell_sites.add_site(LatticeSite(0.5, 0.5, 0.5))
    site_editor.show()

    app.exec_()