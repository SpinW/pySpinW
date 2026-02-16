import json
import os
from enum import Enum
from dataclasses import dataclass, asdict

from PySide6.QtCore import Signal, Qt, QSize
from PySide6.QtGui import QIcon, QPainter
from PySide6.QtWidgets import QWidget, QPushButton, QHBoxLayout, QSpacerItem, QSizePolicy, QLabel, QSlider, QComboBox

from pyspinw.gui.icons.iconload import load_icon, load_pixmap


class RenderMode(Enum):
    default = "default"

@dataclass
class DisplayOptions:
    show_sites: bool = True
    show_couplings: bool = True
    show_anisotropies: bool = True
    show_unit_cell: bool = True
    show_supercell: bool = True

    show_nonmagnetic_atoms: bool = True
    use_atomic_radii: bool = True
    show_atoms_not_moments: bool = False

    atom_moment_scaling: float = 1.0
    coupling_scaling: float = 1.0

    render_mode: RenderMode = RenderMode.default

    def serialise(self) -> str:
        out = asdict(self)
        out["render_mode"] = self.render_mode.value

        return json.dumps(out)

    @staticmethod
    def deserialise(serialised: str):
        data = json.loads(serialised)

        data["render_mode"] = RenderMode(data["render_mode"])

        return DisplayOptions(**data)




class IconWidget(QWidget):
    """ Display an icon in a widget """

    def __init__(self, icon_name: str, alt_text="", parent=None):
        super().__init__(parent)
        self.icon = load_icon(icon_name)
        self.setToolTip(alt_text)
        self.setMinimumSize(QSize(24,24))

    def paintEvent(self, event):
        """ Qt override for painting """
        painter = QPainter(self)
        self.icon.paint(
            painter,
            self.rect(),
            Qt.AlignCenter,
            QIcon.Normal,
            QIcon.Off
        )

class GraphicalSlider(QWidget):
    """ Slider class """

    value_changed = Signal()

    def __init__(self,
                 min_value, max_value, start_value,
                 left_label: QWidget, right_label: QWidget,
                 alt_text: str, parent=None):

        super().__init__(parent)

        self.min_value = min_value
        self.max_value = max_value
        start_position = int(100 * (start_value - min_value) / (max_value - min_value))

        if start_position < 0:
            start_position = 0

        elif start_position > 100:
            start_position = 100

        self.slider = QSlider(Qt.Horizontal)
        self.slider.setMinimum(0)
        self.slider.setMaximum(100)
        self.slider.setValue(start_position)
        self.slider.valueChanged.connect(self.on_slider_changed)
        self.slider.setMinimumSize(QSize(80, 20))

        self.main_layout = QHBoxLayout()

        self.setToolTip(alt_text)
        left_label.setToolTip(alt_text)
        self.slider.setToolTip(alt_text)
        right_label.setToolTip(alt_text)

        self.main_layout.addWidget(left_label)
        self.main_layout.addWidget(self.slider)
        self.main_layout.addWidget(right_label)

        self.setLayout(self.main_layout)
        
    def on_slider_changed(self):
        """ Callback to send signal on """
        self.value_changed.emit()

    def value(self) -> float:
        """ Get the value in floats """
        return self.min_value + 0.01 * self.slider.value() * (self.max_value - self.min_value)

class DisplayOptionsToolbar(QWidget):
    """ Toolbar for setting values """

    settings_filename = "render_settings.conf"

    renderOptionsChanged = Signal()

    def _add_slider(self,
                    min_value: float, max_value: float, start_value: float,
                    left_label: QWidget, right_label: QWidget, alt_text: str):
        """ Add a slider and wire it up """

        slider = GraphicalSlider(
            min_value, max_value, start_value,
            left_label, right_label, alt_text)

        slider.value_changed.connect(self._on_change)
        self.bar_layout.addWidget(slider)

        return slider

    def _add_toggle_button(self, alt_text: str, icon: str | None = None, value: bool = True):
        """ Add a toggle button and wire it up"""

        btn = QPushButton()

        if icon is not None:
            icon = load_icon(icon)
            btn.setIcon(icon)

        btn.setIconSize(QSize(32, 32))
        btn.setMinimumSize(QSize(40, 40))
        btn.setCheckable(True)
        btn.setToolTip(alt_text)
        btn.setChecked(value)

        self.bar_layout.addWidget(btn)
        btn.clicked.connect(self._on_change)

        return btn

    def __init__(self, parent=None):
        super().__init__(parent)

        self.bar_layout = QHBoxLayout()

        #
        # Load settings
        #

        if os.path.exists(self.settings_filename):
            with open(self.settings_filename, 'r') as file:
                settings = DisplayOptions.deserialise(file.read())

        else:
            settings = DisplayOptions()

        #
        # Show hide options
        #

        self.show_sites = self._add_toggle_button(alt_text="Show Sites",
                                                  icon="show_moments",
                                                  value=settings.show_sites)


        self.show_couplings = self._add_toggle_button(alt_text="Show Exchanges",
                                                      icon="show_exchanges",
                                                      value=settings.show_couplings)

        self.show_anisotropies = self._add_toggle_button(alt_text="Show Anisotropies",
                                                         icon="anisotropies",
                                                         value=settings.show_anisotropies)

        self.show_cell = self._add_toggle_button("Show Unit Cell",
                                                 icon="unitcell",
                                                 value=settings.show_unit_cell)

        self.show_supercell = self._add_toggle_button("Show Supercell",
                                                      icon="supercell",
                                                      value=settings.show_supercell)


        self.bar_layout.addSpacerItem(QSpacerItem(20, 0, QSizePolicy.Minimum, QSizePolicy.Minimum))

        #
        # Various visual properties
        #

        self.atom_or_moments = self._add_toggle_button("Switch between showing atoms and moments",
                                                       icon="momentatoms",
                                                       value=settings.show_atoms_not_moments)


        self.show_nonmagnetic = self._add_toggle_button("Show non-magnetic sites",
                                                        icon="nonmagnetic",
                                                        value=settings.show_nonmagnetic_atoms)


        self.moment_scale_slider = self._add_slider(
            0, 1, settings.atom_moment_scaling,
            left_label=IconWidget("small_moments", "Smaller Sites"),
            right_label=IconWidget("big_moments", "Larger Sites"),
            alt_text="Site scale factor")

        self.scale_atoms = self._add_toggle_button("Scale atoms by atomic radii",
                                                   icon="atomsizes",
                                                   value=settings.use_atomic_radii)

        self.coupling_scale_slider = self._add_slider(
            0, 1, settings.coupling_scaling,
            left_label=IconWidget("small_exchange", "Smaller Exchanges"),
            right_label=IconWidget("big_exchange", "Larger Exchanges"),
            alt_text="Exchange thickness")



        self.bar_layout.addSpacerItem(QSpacerItem(20, 0, QSizePolicy.Minimum, QSizePolicy.Minimum))

        #
        # Render mode
        #

        self.render_mode_combo = QComboBox()
        for item in RenderMode:
            self.render_mode_combo.addItem(item.value)

        self.render_mode_combo.setCurrentText(settings.render_mode.value)

        self.render_mode_combo.setAccessibleDescription("Rendering method")
        self.bar_layout.addWidget(self.render_mode_combo)


        # Pad right, and set layout
        self.bar_layout.addSpacerItem(QSpacerItem(0, 0, QSizePolicy.Expanding, QSizePolicy.Minimum))
        self.setLayout(self.bar_layout)

    def _on_change(self):
        """ Called when anything changes, send signal """
        self.renderOptionsChanged.emit()

    def display_options(self) -> DisplayOptions:
        """ Get the current render options """

        render_mode = RenderMode(self.render_mode_combo.currentText())

        return DisplayOptions(
            show_sites = self.show_sites.isChecked(),
            show_couplings = self.show_couplings.isChecked(),
            show_anisotropies = self.show_anisotropies.isChecked(),
            show_unit_cell = self.show_cell.isChecked(),
            show_supercell = self.show_supercell.isChecked(),
            show_nonmagnetic_atoms = self.show_nonmagnetic.isChecked(),
            use_atomic_radii = self.scale_atoms.isChecked(),
            atom_moment_scaling = self.moment_scale_slider.value(),
            coupling_scaling = self.coupling_scale_slider.value(),
            render_mode = render_mode)

    def save_settings(self):
        with open(self.settings_filename, 'w') as file:
            file.write(self.display_options().serialise())

    def closeEvent(self, event):
        self.save_settings()

        super().closeEvent(event)

if __name__ == "__main__":
    import sys
    from PySide6.QtWidgets import QApplication

    app = QApplication()
    widget = DisplayOptionsToolbar()
    widget.resize(800, 600)
    widget.show()
    sys.exit(app.exec())
