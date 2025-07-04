""" Miscellaneous widgets that are generally useful"""

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import QLabel, QWidget, QSlider, QVBoxLayout


class QRightLabel(QLabel):
    """ Right aligned label"""

    def __init__(self, text: str):
        super().__init__(text)
        self.setAlignment(Qt.AlignRight | Qt.AlignVCenter)

class QLeftLabel(QLabel):
    """ Right aligned label"""

    def __init__(self, text: str):
        super().__init__(text)
        self.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)



class FloatSlider(QWidget):
    """ Floating point slider """

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
        """ Slider value as a float"""
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


