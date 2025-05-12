from dataclasses import dataclass

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QVBoxLayout, QWidget, QComboBox

from pyspinw.symmetry.system import crystal_systems, crystal_system_name_lookup


class SymmetryWidget(QWidget):

    symmetry_changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent=parent)

        layout = QVBoxLayout()

        self.system_combo = QComboBox(self)
        self.system_combo.currentTextChanged.connect(self._on_crystal_system_changed)

        for crystal_system in crystal_systems:
            self.system_combo.addItem(crystal_system.name)

        self.system_combo.setCurrentText("Triclinic")

        layout.addWidget(self.system_combo)
        self.setLayout(layout)



    def _on_crystal_system_changed(self):
        self.symmetry_changed.emit()

    @property
    def current_crystal_system(self):
        return crystal_system_name_lookup[self.system_combo.currentText()]

