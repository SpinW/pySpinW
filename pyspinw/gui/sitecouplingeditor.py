from PySide6.QtCore import Signal
from PySide6.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QSpacerItem, QSizePolicy

from pyspinw.gui.couplingtable import CouplingTable
from pyspinw.gui.couplingcreator import CouplingCreatorWindow
from pyspinw.gui.crystalviewer.actionlabel import ActionLabel
from pyspinw.gui.helperwidgets.dockwidget import SpinWDockWidget
from pyspinw.gui.sitetable import SiteTable
from pyspinw.gui.symmetry_settings import SymmetrySettings, DEFAULT_SYMMETRY
from pyspinw.site import LatticeSite


class SiteButtons(QWidget):

    add = Signal()
    remove = Signal()
    reify = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)

        layout = QHBoxLayout()

        self.add_button = QPushButton("Add")
        self.remove_button = QPushButton("Remove")
        self.reify_button = QPushButton("Make Inequivalent")

        self.add_button.clicked.connect(self._on_add)
        self.remove_button.clicked.connect(self._on_remove)
        self.reify_button.clicked.connect(self._on_reify)

        self.remove_button.setEnabled(False)
        self.reify_button.setEnabled(False)

        layout.addWidget(self.add_button)
        layout.addWidget(self.remove_button)
        layout.addWidget(self.reify_button)

        self.setLayout(layout)


    def _on_add(self):
        self.add.emit()

    def _on_remove(self):
        self.remove.emit()

    def _on_reify(self):
        self.reify.emit()

class CouplingButtons(QWidget):

    add_clicked = Signal()
    remove_clicked = Signal()

    def __init__(self, parent=None):
        super().__init__(parent=parent)

        self.add_button = QPushButton("Add")
        self.remove_button = QPushButton("Remove")

        self.add_button.clicked.connect(self._on_add)
        self.remove_button.clicked.connect(self._on_remove)

        layout = QHBoxLayout()

        layout.addWidget(self.add_button)
        layout.addWidget(self.remove_button)

        self.setLayout(layout)

    def _on_add(self):
        """ Add button pressed"""
        self.add_clicked.emit()

    def _on_remove(self):
        """ Remove button pressed"""
        self.remove_clicked.emit()

class SiteAndCouplingEditor(SpinWDockWidget):
    """ Editor dock window for magnetic sites"""

    graphics_relevant_change = Signal()

    def __init__(self, parent=None):
        super().__init__(parent=parent)

        self._symmetry = DEFAULT_SYMMETRY


        self.setWindowTitle("Sites")

        widget = QWidget()
        layout = QVBoxLayout()

        self.setWidget(widget)

        self.coupling_window = CouplingCreatorWindow()

        self.site_table = SiteTable(symmetry=self._symmetry)
        self.fix_bad_sites_label = ActionLabel(
            "The moments on some sites appear to conflict, this is either because "
            "your symmetry is too high for your system, or because the moments should"
            "be zero.",
            action_text="Click to zero bad moments.")

        self.site_buttons = SiteButtons()

        self.coupling_table = CouplingTable()
        self.coupling_buttons = CouplingButtons()

        self.coupling_details_container = QWidget()

        widget.setLayout(layout)

        layout.addWidget(QLabel("Sites"))

        layout.addWidget(self.site_table)
        layout.addWidget(self.fix_bad_sites_label)
        layout.addWidget(self.site_buttons)

        layout.addWidget(QLabel("Couplings"))

        layout.addWidget(self.coupling_table)
        layout.addWidget(self.coupling_buttons)
        layout.addWidget(self.coupling_details_container)

        layout.addSpacerItem(QSpacerItem(10, 10, QSizePolicy.Minimum, QSizePolicy.Expanding))

        # Site connections

        self.site_buttons.add.connect(self._on_add)

        self.site_table.graphics_relevant_change.connect(self._on_graphics_relevant_change)
        self.site_table.magnetic_symmetry_broken.connect(self._on_magnetic_symmetry_broken)
        self.site_table.magnetic_symmetry_ok.connect(self._on_magnetic_symmetry_ok)

        self.site_table.implied_selected_state_changed.connect(self._on_implied_selected_state_changed)
        self.site_table.non_implied_selected_state_changed.connect(self._on_non_implied_selected_state_changed)

        self.fix_bad_sites_label.action.connect(self.site_table.magnetic_symmetry_autofix)

        # Coupling connections
        self.coupling_buttons.add_clicked.connect(self._on_add_coupling_clicked)
        self.coupling_window.couplings_accepted.connect(self._on_new_couplings_accepted)

    @property
    def symmetry(self):
        return self._symmetry

    @symmetry.setter
    def symmetry(self, symmetry: SymmetrySettings):
        self._symmetry = symmetry
        self.site_table.symmetry = symmetry

        self.update()

    def _on_add(self):
        self.site_table.add_site(LatticeSite.create(0,0,0,0,0,1,name="New Site"))

    def _on_graphics_relevant_change(self):
        self.graphics_relevant_change.emit()

    def _on_magnetic_symmetry_broken(self):
        self.fix_bad_sites_label.setVisible(True)

    def _on_magnetic_symmetry_ok(self):
        self.fix_bad_sites_label.setVisible(False)

    def _on_non_implied_selected_state_changed(self, state: bool):
        self.site_buttons.remove_button.setEnabled(state)

    def _on_implied_selected_state_changed(self, state: bool):
        self.site_buttons.reify_button.setEnabled(state)

    @property
    def sites_for_drawing(self):
        return self.site_table.sites_for_drawing

    def _on_add_coupling_clicked(self):
        self.coupling_window.show_creator(self.site_table.selected_sites, self.symmetry.unit_cell)

    def _on_new_couplings_accepted(self):
        self.coupling_table.add_couplings(self.coupling_window.couplings)
