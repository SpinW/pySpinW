from PySide6.QtGui import QDoubleValidator
from PySide6.QtWidgets import QLineEdit, QWidget, QHBoxLayout, QApplication
from PySide6.QtCore import Signal

class FloatField(QWidget):

    updated = Signal(float)

    def __init__(self, value: float):
        super().__init__()

        self.text_field = QLineEdit(str(value))

        validator = QDoubleValidator(0.0, 100.0, 2, self)
        validator.setNotation(QDoubleValidator.StandardNotation)
        self.text_field.setValidator(validator)

        layout = QHBoxLayout()
        layout.addWidget(self.text_field)

        self.setLayout(layout)

if __name__ == "__main__":

    app = QApplication([])

    field = FloatField(17.84763)

    field.show()

    app.exec_()
