from PySide6.QtCore import Signal, Qt
from PySide6.QtWidgets import QGridLayout, QWidget, QComboBox, QLabel
from ase.lattice import bravais_lattices

from pyspinw.symmetry.system import lattice_systems, lattice_system_name_lookup
from pyspinw.symmetry.bravais import lattice_type_name_lookup
from pyspinw.symmetry.group import spacegroup_lattice_symbol_lookup, SymmetryGroup, MagneticSpaceGroup, SpaceGroup, \
    spacegroup_symbol_lookup

from pyspinw.gui.helperwidgets.misc import QRightLabel

class SymmetryWidget(QWidget):

    symmetry_changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent=parent)

        layout = QGridLayout()

        self.system_combo = QComboBox(self)
        self.bravais_type_combo = QComboBox(self)
        self.space_group_combo = QComboBox(self)
        self.magnetic_spacegroup_combo = QComboBox(self)


        # Fill in the crystal system combo, and set
        for crystal_system in lattice_systems:
            self.system_combo.addItem(crystal_system.name)

        self.system_combo.setCurrentText("Triclinic")

        # Fill in the bravais lattice combo
        self._set_bravais_combo()
        self._set_spacegroup_combo()
        self._set_magnetic_spacegroup_combo()

        self.system_combo.currentTextChanged.connect(self._on_lattice_system_changed)
        self.bravais_type_combo.currentTextChanged.connect(self._on_bravais_changed)
        self.space_group_combo.currentTextChanged.connect(self._on_spacegroup_changed)
        self.magnetic_spacegroup_combo.currentTextChanged.connect(self._on_magnetic_spacegroup_changed)

        # Do the layout

        layout.addWidget(QRightLabel("Lattice System"), 0, 0)
        layout.addWidget(self.system_combo, 0, 1)

        layout.addWidget(QRightLabel("Bravais Type"), 1, 0)
        layout.addWidget(self.bravais_type_combo, 1, 1)

        layout.addWidget(QRightLabel("Spacegroup"), 2, 0)
        layout.addWidget(self.space_group_combo, 2, 1)

        layout.addWidget(QRightLabel("Magnetic Group"), 3, 0)
        layout.addWidget(self.magnetic_spacegroup_combo, 3, 1)

        self.setLayout(layout)

    def _set_bravais_combo(self):
        current_selection = self.bravais_type_combo.currentText()

        self.bravais_type_combo.blockSignals(True)
        self.bravais_type_combo.clear()

        names = [bravais.name for bravais in self.current_lattice_system.bravais_options.bravias]

        for name in names:
            self.bravais_type_combo.addItem(name)

        if current_selection in names:
            self.bravais_type_combo.setCurrentText(current_selection)
        else:
            self.bravais_type_combo.setCurrentText(names[0])

        self.bravais_type_combo.blockSignals(False)

    def _set_spacegroup_combo(self):
        current_selection = self.space_group_combo.currentText()

        self.space_group_combo.blockSignals(True)

        symbol = self.current_bravais_symbol
        names = [spacegroup.symbol for spacegroup in spacegroup_lattice_symbol_lookup[symbol]]

        self.space_group_combo.clear()
        for name in names:
            self.space_group_combo.addItem(name)

        if current_selection in names:
            self.space_group_combo.setCurrentText(current_selection)
        else:
            self.space_group_combo.setCurrentText(names[0])

        self.space_group_combo.blockSignals(False)

    def _set_magnetic_spacegroup_combo(self):
        current_selection = self.magnetic_spacegroup_combo.currentText()

        names = [magnetic.symbol for magnetic in self.current_spacegroup.magnetic_variants]

        self.magnetic_spacegroup_combo.clear()
        for name in names:
            self.magnetic_spacegroup_combo.addItem(name)

        if current_selection in names:
            self.magnetic_spacegroup_combo.setCurrentText(current_selection)
        else:
            self.magnetic_spacegroup_combo.setCurrentText(names[0])


    def _on_lattice_system_changed(self):
        self._set_bravais_combo()
        self._set_spacegroup_combo()
        self._set_magnetic_spacegroup_combo()

        # Do this last
        self.symmetry_changed.emit()

    def _on_bravais_changed(self):
        self._set_spacegroup_combo()
        self._set_magnetic_spacegroup_combo()

        # Do this last
        self.symmetry_changed.emit()

    def _on_spacegroup_changed(self):
        self._set_magnetic_spacegroup_combo()

        # Do this last
        self.symmetry_changed.emit()

    def _on_magnetic_spacegroup_changed(self):

        # Do this last (not that theres currently anything before)
        self.symmetry_changed.emit()

    @property
    def current_lattice_system(self):
        return lattice_system_name_lookup[self.system_combo.currentText()]

    @property
    def current_bravais(self):
        return lattice_type_name_lookup[self.bravais_type_combo.currentText()]


    @property
    def current_bravais_symbol(self):
        lattice = self.current_lattice_system
        bravais = self.current_bravais
        return lattice.letter + bravais.letter

    @property
    def current_spacegroup_name(self):
        return self.space_group_combo.currentText()

    @property
    def current_spacegroup(self) -> SpaceGroup:
        return spacegroup_symbol_lookup[self.current_spacegroup_name]

    def lattice_autoset(self, lattice_system_name):
        if lattice_system_name in lattice_system_name_lookup: # Check for existence
            self.system_combo.setCurrentText(lattice_system_name)
        else:
            raise ValueError("Not a valid LatticeSystem name")