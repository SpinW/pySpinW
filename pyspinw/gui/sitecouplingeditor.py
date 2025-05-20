from PySide6.QtCore import Signal
from PySide6.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QSpacerItem, QSizePolicy

from pyspinw.gui.decorated import DecoratedSite
from pyspinw.gui.helperwidgets.dockwidget import SpinWDockWidget
from pyspinw.gui.helperwidgets.alternatesitetable import SiteTable
from pyspinw.gui.symmetry_settings import SymmetrySettings
from pyspinw.site import LatticeSite


class SiteButtons(QWidget):

    add = Signal()
    remove = Signal()
    reify_button = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)

        layout = QHBoxLayout()

        self.add_button = QPushButton("Add")
        self.remove_button = QPushButton("Remove")
        self.reify_button = QPushButton("Make Explicit")

        self.add_button.clicked.connect(self._on_add)
        self.remove_button.clicked.connect(self._on_remove)
        self.reify_button.clicked.connect(self._on_reify)

        layout.addWidget(self.add_button)
        layout.addWidget(self.remove_button)
        layout.addWidget(self.reify_button)

        self.setLayout(layout)


    def _on_add(self):
        self.add.emit()

    def _on_remove(self):
        self.remove.emit()

    def _on_reify(self):
        self.reify_button.emit()

class SiteAndCouplingEditor(SpinWDockWidget):
    """ Editor dock window for magnetic sites"""

    graphics_relevant_change = Signal()

    def __init__(self, parent=None):
        super().__init__(parent=parent)

        self.setWindowTitle("Sites")

        widget = QWidget()
        layout = QVBoxLayout()

        self.setWidget(widget)

        self.site_table = SiteTable()
        self.site_buttons = SiteButtons()

        widget.setLayout(layout)

        layout.addWidget(QLabel("Sites"))
        layout.addWidget(self.site_table)
        layout.addWidget(self.site_buttons)

        layout.addSpacerItem(QSpacerItem(10, 10, QSizePolicy.Minimum, QSizePolicy.Expanding))

        self.site_buttons.add.connect(self._on_add)

        self._symmetry: None

        self.site_table.graphics_relevant_change.connect(self._on_graphics_relevant_change)

    @property
    def symmetry(self):
        return self._symmetry

    @symmetry.setter
    def symmetry(self, symmetry: SymmetrySettings):
        self._symmetry = symmetry
        self.site_table.symmetry = symmetry

        self.update()

    def _on_add(self):
        self.site_table.add_site(LatticeSite(0,0,0,0,0,1,name="New Site"))

    def _on_graphics_relevant_change(self):
        self.graphics_relevant_change.emit()

    @property
    def sites_for_drawing(self):
        return self.site_table.sites_for_drawing