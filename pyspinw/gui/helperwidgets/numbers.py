from PySide6.QtGui import QDoubleValidator
from PySide6.QtWidgets import QLineEdit, QWidget, QHBoxLayout, QApplication
from PySide6.QtCore import Signal

import sys

class FloatField(QWidget):


    changed = Signal(float)

    def __init__(self,
                 value: float,
                 bottom: float | None = None,
                 top: float | None = None,
                 show_slider: bool = True):
        """Floating point field wrapper

        :param value: initial value
        :param bottom: lower bound
        :param top: upper bound
        :param show_slider: turn on the slider
        """
        super().__init__()

        self.text_field = QLineEdit(str(value))

        self.bottom = bottom
        self.top = top

        if bottom is None:
            bottom = -sys.float_info.max

        if top is None:
            top = sys.float_info.max

        validator = QDoubleValidator(self, bottom=bottom, top=top, decimals=100)
        validator.setNotation(QDoubleValidator.StandardNotation)
        self.text_field.setValidator(validator)

        layout = QHBoxLayout()
        layout.addWidget(self.text_field)

        self.text_field.textChanged.connect(self._on_changed)

        self.setLayout(layout)

        # Remove padding
        self.setContentsMargins(0,0,0,0)


        self._value = value

    def _on_changed(self):
        """Event when text is changed"""

        try:
            potential_value = float(self.text_field.text())

            if self.bottom is None or potential_value >= self.bottom and \
                    self.top is None or potential_value <= self.top:

                self._value = potential_value

        except ValueError:
            pass

        self.changed.emit(self.value)

    @property
    def value(self) -> float:
        """Floating point value of field"""
        return self._value

if __name__ == "__main__":

    app = QApplication([])

    field = FloatField(17.84763)

    field.show()

    app.exec_()
