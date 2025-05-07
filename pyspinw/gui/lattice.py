from PySide6.QtWidgets import QBoxLayout, QGridLayout, QWidget, QApplication, QLabel

from pyspinw.gui.helperwidgets.numbers import FloatField
from pyspinw.unitcell import UnitCell


class LatticeParameters(QWidget):
    def __init__(self, unit_cell: UnitCell | None = None, parent=None):
        super().__init__(parent)

        self.a = FloatField(1.0, bottom=0.0)
        self.b = FloatField(1.0, bottom=0.0)
        self.c = FloatField(1.0, bottom=0.0)

        self.alpha = FloatField(90, 0, 180)
        self.beta = FloatField(90, 0, 180)
        self.gamma = FloatField(90, 0, 180)

        layout = QGridLayout(parent=self)
        # layout.setVerticalSpacing(0)
        self.setLayout(layout)

        layout.addWidget(QLabel("a"), 0, 0)
        layout.addWidget(QLabel("b"), 1, 0)
        layout.addWidget(QLabel("c"), 2, 0)

        layout.addWidget(QLabel("α"), 3, 0)
        layout.addWidget(QLabel("β"), 4, 0)
        layout.addWidget(QLabel("γ"), 5, 0)

        layout.addWidget(self.a, 0, 1)
        layout.addWidget(self.b, 1, 1)
        layout.addWidget(self.c, 2, 1)

        layout.addWidget(self.alpha, 3, 1)
        layout.addWidget(self.beta, 4, 1)
        layout.addWidget(self.gamma, 5, 1)

        layout.addWidget(QLabel("Å"), 0, 2)
        layout.addWidget(QLabel("Å"), 1, 2)
        layout.addWidget(QLabel("Å"), 2, 2)
        layout.addWidget(QLabel("°"), 3, 2)
        layout.addWidget(QLabel("°"), 4, 2)
        layout.addWidget(QLabel("°"), 5, 2)


if __name__ == "__main__":
    app = QApplication([])

    parameter_widget = LatticeParameters()
    parameter_widget.show()

    app.exec_()
