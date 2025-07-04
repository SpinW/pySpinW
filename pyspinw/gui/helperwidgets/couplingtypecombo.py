""" Combo box for coupling"""

from PySide6.QtWidgets import QComboBox

from pyspinw.coupling import couplings


class CouplingTypeCombo(QComboBox):
    """ Combo box for coupling types """

    def __init__(self, coupling_type: str, parent=None):
        super().__init__(parent)

        for coupling in couplings:
            self.addItem(coupling.coupling_type)

        self.setCurrentText(coupling_type)
