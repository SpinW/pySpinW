from PySide6.QtWidgets import QWidget, QVBoxLayout, QTreeWidget, QLabel, QApplication, QTreeWidgetItem, \
    QAbstractItemView, QListWidget, QListWidgetItem, QHBoxLayout, QSpacerItem, QSizePolicy, QPushButton

from pyspinw.gui.coordinateeditor import CoordinateEditor
from pyspinw.gui.helperwidgets.dockwidget import SpinWDockWidget
from pyspinw.gui.helperwidgets.sitetable import SiteTable
from pyspinw.site import LatticeSite
from pyspinw.symmetry.group import SymmetryGroup


class AsymmetricSiteItem(QListWidgetItem):
    def __init__(self, site: LatticeSite):
        super().__init__("A site")

        self.site = site


class AsymmetricButtons(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent=parent)

        layout = QHBoxLayout()
        layout.addSpacerItem(QSpacerItem(10, 0,  QSizePolicy.Expanding, QSizePolicy.Minimum))

        self.add_button = QPushButton("Add")
        self.delete_button = QPushButton("Remove")

        layout.addWidget(self.add_button)
        layout.addWidget(self.delete_button)

        self.setLayout(layout)


class ImpliedButtons(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent=parent)

        layout = QHBoxLayout()
        layout.addSpacerItem(QSpacerItem(10, 0,  QSizePolicy.Expanding, QSizePolicy.Minimum))

        self.explicit = QPushButton("Make Explicit")

        layout.addWidget(self.explicit)

        self.setLayout(layout)


class SiteEditor(SpinWDockWidget):
    """ Editor dock window for magnetic sites"""
    def __init__(self, parent=None):
        super().__init__(parent=parent)

        self.setWindowTitle("Sites")

        widget = QWidget()
        layout = QVBoxLayout()

        self.setWidget(widget)

        self.asymmetric_cell_sites = SiteTable()
        self.asymmetric_cell_buttons = AsymmetricButtons()
        self.implied_sites = SiteTable()
        self.implied_buttons = ImpliedButtons()
        self.coordinate_editor = CoordinateEditor()

        widget.setLayout(layout)

        layout.addWidget(QLabel("Asymmetric Cell"))
        layout.addWidget(self.asymmetric_cell_sites)
        layout.addWidget(self.asymmetric_cell_buttons)

        layout.addWidget(QLabel("Implied Sites"))
        layout.addWidget(self.implied_sites)
        layout.addWidget(self.implied_buttons)

        layout.addWidget(QLabel("Values"))
        layout.addWidget(self.coordinate_editor)

        layout.addSpacerItem(QSpacerItem(10, 10, QSizePolicy.Minimum, QSizePolicy.Expanding))

    # def update(self, asymmetric_sites: list[LatticeSite], symmetry: SymmetryGroup):
    #     """ Update the implicit sites given by the symmetry """
    #
    #     # Does this update affect what sites there are, rather than the locations, if so, we want to give a warning
    #     # We don't want to warn if we've deleted/added something from the initial list, just if its moved to a special point
    #     # Or if we have changed the symmetry
    #     # We should also warn if we're deleting something that affects a constraint



if __name__ == "__main__":
    app = QApplication([])

    site_editor = SiteEditor()
    site_editor.asymmetric_cell_sites.add_site(LatticeSite(0.5, 0.5, 0.5))
    site_editor.asymmetric_cell_sites.add_site(LatticeSite(0.5, 0.5, 0.5))
    site_editor.asymmetric_cell_sites.add_site(LatticeSite(0.5, 0.5, 0.5))
    site_editor.show()

    app.exec_()