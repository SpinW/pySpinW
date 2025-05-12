from PySide6.QtCore import Signal, Qt
from PySide6.QtWidgets import QGridLayout, QWidget, QComboBox, QLabel

from pyspinw.symmetry.system import crystal_systems, crystal_system_name_lookup

class QRightLabel(QLabel):
    def __init__(self, text: str):
        super().__init__(text)
        self.setAlignment(Qt.AlignRight)

class SymmetryWidget(QWidget):

    symmetry_changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent=parent)

        layout = QGridLayout()

        self.system_combo = QComboBox(self)
        self.bravais_type_combo = QComboBox(self)
        self.space_group_combo = QComboBox(self)
        self.magnetic_space_group_combo = QComboBox(self)

        self.system_combo.currentTextChanged.connect(self._on_crystal_system_changed)
        self.bravais_type_combo.currentTextChanged.connect(self._on_bravais_changed)

        # Fill in the crystal system combo, and set
        for crystal_system in crystal_systems:
            self.system_combo.addItem(crystal_system.name)

        self.system_combo.setCurrentText("Triclinic")

        # Fill in the bravais lattice combo
        self._set_bravais_combo()

        # Do the layout

        layout.addWidget(QRightLabel("Crystal System"), 0, 0)
        layout.addWidget(self.system_combo, 0, 1)

        layout.addWidget(QRightLabel("Bravais Type"), 1, 0)
        layout.addWidget(self.bravais_type_combo, 1, 1)

        layout.addWidget(QRightLabel("Spacegroup"), 2, 0)
        layout.addWidget(self.space_group_combo, 2, 1)

        layout.addWidget(QRightLabel("Magnetic Group"), 3, 0)
        layout.addWidget(self.magnetic_space_group_combo, 3, 1)

        self.setLayout(layout)

    def _set_bravais_combo(self):
        current_selection = self.bravais_type_combo.currentText()

        self.bravais_type_combo.clear()

        names = [bravais.name for bravais in self.current_crystal_system.bravais_options.bravias]

        for name in names:
            self.bravais_type_combo.addItem(name)

        if current_selection in names:
            self.bravais_type_combo.setCurrentText(current_selection)

    def _on_crystal_system_changed(self):
        self._set_bravais_combo()
        
        # Do this last
        self.symmetry_changed.emit()

    def _on_bravais_changed(self):
        pass

    @property
    def current_crystal_system(self):
        return crystal_system_name_lookup[self.system_combo.currentText()]

