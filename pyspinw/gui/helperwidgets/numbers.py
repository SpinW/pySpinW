from PySide6.QtGui import QDoubleValidator
from PySide6.QtWidgets import QLineEdit, QWidget, QHBoxLayout, QApplication, QSlider, QPushButton, QVBoxLayout
from PySide6.QtCore import Signal, Qt

import sys

class FloatSlider(QWidget):

    changed = Signal(float)
    resolution = 1000

    def __init__(self, orientation: Qt.Orientation, value: float, bottom: float, top: float, parent=None):
        super().__init__(parent)

        self.slider = QSlider(orientation)
        self.slider.setRange(0, FloatSlider.resolution)

        self.top = top
        self.bottom = bottom
        self.scale = (top - bottom) / FloatSlider.resolution

        self._value = value

        layout = QVBoxLayout()
        layout.addWidget(self.slider)
        self.setLayout(layout)

        self._update_slider(value)


        self.slider.sliderMoved.connect(self._on_changed)

    def _on_changed(self):
        self.changed.emit(self.value)

    @property
    def value(self) -> float:
        return self.slider.value() * self.scale + self.bottom

    @value.setter
    def value(self, value):
        self._value = value
        self._update_slider(value, no_message=True)

    def _update_slider(self, value: float, no_message: bool=False):
        slider_value = max(min(int((value - self.bottom) / self.scale), FloatSlider.resolution), 0)

        if no_message:
            self.slider.blockSignals(True)

        self.slider.setValue(slider_value)

        if no_message:
            self.slider.blockSignals(False)




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
        self.slider_toggle = QPushButton("...")

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

        outer_layout = QVBoxLayout()
        inner_widget = QWidget(self)

        inner_layout = QHBoxLayout()
        inner_widget.setLayout(inner_layout)
        inner_layout.addWidget(self.text_field)

        outer_layout.addWidget(inner_widget)

        if allow_slider:
            inner_layout.addWidget(self.slider_toggle)
            outer_layout.addWidget(self.slider)

        self.setLayout(outer_layout)

        # Remove padding
        # self.setContentsMargins(0,0,0,0)

        self._value = value

        # self.setStyleSheet("background-color: lightblue;")

    def _on_slider_changed(self):
        """ Event when the slider is changed"""

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

if __name__ == "__main__":

    app = QApplication([])

    field = FloatField(17.84763, slider_bottom=0, slider_top=100)

    field.show()

    app.exec_()
