import sys

from PySide6.QtGui import QDoubleValidator
from PySide6.QtWidgets import QLineEdit, QWidget, QHBoxLayout, QApplication, QSlider, QPushButton, QVBoxLayout
from PySide6.QtCore import Signal, Qt

from pyspinw.gui.helperwidgets.misc import FloatSlider


class FloatField(QWidget):

    changed = Signal(float)

    def __init__(self,
                 value: float,
                 bottom: float | None = None,
                 top: float | None = None,
                 allow_slider: bool = True,
                 minimise_slider: bool = True,
                 slider_bottom: float | None = None,
                 slider_top: float | None = None,
                 parent=None):

        """Floating point field wrapper

        :param value: initial value
        :param bottom: lower bound
        :param top: upper bound
        :param allow_slider: turn on the slider
        :param minimise_slider: make the slider small
        """
        super().__init__(parent=parent)


        self.bottom = bottom
        self.top = top
        self._allow_slider = allow_slider

        if allow_slider:

            if slider_bottom is None:
                slider_bottom = bottom

            if slider_top is None:
                slider_top = top

            if slider_bottom is None:
                raise ValueError("No bottom of range information for slider")

            if slider_top is None:
                raise ValueError("No top of range information for slider")

            self.slider = FloatSlider(Qt.Orientation.Horizontal, value, slider_bottom, slider_top)
            self.slider.changed.connect(self._on_slider_changed)

        else:
            self.slider = QSlider(Qt.Orientation.Horizontal)

        if minimise_slider:
            self.slider.setVisible(False)

        self.text_field = QLineEdit(str(value))
        self.slider_toggle = QPushButton("◄ ►", flat=True)
        self.slider_toggle.setMaximumWidth(30)

        if bottom is None:
            bottom = -sys.float_info.max

        if top is None:
            top = sys.float_info.max

        validator = QDoubleValidator(self, bottom=bottom, top=top, decimals=100)
        validator.setNotation(QDoubleValidator.StandardNotation)
        self.text_field.setValidator(validator)

        # Slots

        self.text_field.textChanged.connect(self._on_text_changed)
        self.slider_toggle.clicked.connect(self.toggle_slider)

        # Layout


        inner_widget = QWidget(parent=self)
        inner_layout = QHBoxLayout()
        inner_layout.setContentsMargins(0, 0, 0, 0)
        inner_layout.setSpacing(0)

        inner_widget.setLayout(inner_layout)
        inner_layout.addWidget(self.text_field)

        outer_layout = QVBoxLayout()
        outer_layout.addWidget(inner_widget)
        outer_layout.setContentsMargins(0, 0, 0, 0)
        outer_layout.setSpacing(0)

        self.setLayout(outer_layout)

        if allow_slider:
            inner_layout.addWidget(self.slider_toggle)
            outer_layout.addWidget(self.slider)

        self._value = value


    def _on_slider_changed(self):
        """ Event when the slider is changed """

        self._value = self.slider.value
        self.text_field.blockSignals(True)
        self.text_field.setText(f"{self._value:.6g}")
        self.text_field.blockSignals(False)

        self.changed.emit(self._value)


    def _on_text_changed(self):
        """Event when text is changed"""

        try:
            potential_value = float(self.text_field.text())

            if self.bottom is None or potential_value >= self.bottom and \
                    self.top is None or potential_value <= self.top:

                self._value = potential_value

        except ValueError:
            pass

        if self._allow_slider:
            self.slider.value = self._value

        self.changed.emit(self.value)

    def toggle_slider(self):
        self.slider.setVisible(not self.slider.isVisible())

    @property
    def value(self) -> float:
        """Floating point value of field"""
        return self._value

    @value.setter
    def value(self, value):
        self._value = value
        self.text_field.setText(str(value))

    def get_value(self):
        """ Value getter that is not a property (needed to avoid certain reference issues)"""
        return self._value

if __name__ == "__main__":

    app = QApplication([])

    field = FloatField(17.84763, slider_bottom=0, slider_top=100, allow_slider=False)


    field.show()

    app.exec_()
