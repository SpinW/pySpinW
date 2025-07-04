""" Index field helper widget (ijk entries) """

from PySide6.QtCore import Signal
from PySide6.QtGui import QIntValidator
from PySide6.QtWidgets import QLineEdit, QWidget, QHBoxLayout


class IndexField(QWidget):
    """ Field for entering 3 indices"""

    value_changed = Signal()

    def __init__(self, value: tuple[int, int, int] = (0,0,0), parent=None):
        super().__init__(parent)

        self.i = QLineEdit(text=str(value[0]), parent=self)
        self.j = QLineEdit(text=str(value[1]), parent=self)
        self.k = QLineEdit(text=str(value[2]), parent=self)

        i_validator = QIntValidator(-1000, 1000, self.i)
        j_validator = QIntValidator(-1000, 1000, self.j)
        k_validator = QIntValidator(-1000, 1000, self.k)

        self.i.setValidator(i_validator)
        self.j.setValidator(j_validator)
        self.k.setValidator(k_validator)

        self.i.textChanged.connect(self._on_value_changed)
        self.j.textChanged.connect(self._on_value_changed)
        self.k.textChanged.connect(self._on_value_changed)

        layout = QHBoxLayout()

        layout.addWidget(self.i)
        layout.addWidget(self.j)
        layout.addWidget(self.k)

        self.setLayout(layout)

    def _on_value_changed(self):
        self.value_changed.emit()

    @property
    def value(self) -> tuple[int, int, int]:
        """ Current value as a tuple of ints"""
        return int(self.i.text()), int(self.j.text()), int(self.k.text())

    @value.setter
    def value(self, value: tuple[int, int, int]):
        self.i.setText(value[0])
        self.j.setText(value[1])
        self.k.setText(value[2])
