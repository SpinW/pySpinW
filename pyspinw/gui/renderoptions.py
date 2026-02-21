""" Toolbar and display options data for the viewer"""
import ctypes
import json
import os

from dataclasses import dataclass, asdict

from PySide6.QtCore import Signal, Qt, QSize
from PySide6.QtGui import QIcon, QPainter
from PySide6.QtWidgets import QWidget, QPushButton, QHBoxLayout, QSpacerItem, QSizePolicy, QSlider

from pyspinw.gui.icons.iconload import load_icon


@dataclass
class DisplayOptions:
    """ Options for how things should display """

    show_sites: bool = True
    show_couplings: bool = True
    show_anisotropies: bool = True
    show_unit_cell: bool = False
    show_supercell: bool = False

    prettify: bool = True

    show_nonmagnetic_atoms: bool = True
    use_atomic_radii: bool = True
    show_atoms_not_moments: bool = False

    atom_moment_scaling: float = 0.35
    coupling_scaling: float = 0.20

    show_cartesian_axes: bool = True
    show_lattice_axes: bool = False
    orthogonal_lattice_axes: bool = False


    def serialise(self) -> str:
        """ Serialise settings object to a string"""
        out = asdict(self)

        return json.dumps(out)

    @staticmethod
    def deserialise(serialised: str):
        """ Deserialise a setting object from a string"""
        try:
            data = json.loads(serialised)

            return DisplayOptions(**data)

        except Exception: # Anything at all goes wrong, just return default

            return DisplayOptions()




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

    displayOptionsChanged = Signal()
    requestViewReset = Signal()

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

    def _toolbar_button(self, alt_text, icon: str | None):
        """ Create a toolbar button with standardised formatting"""
        btn = QPushButton()

        if icon is not None:
            icon = load_icon(icon)
            btn.setIcon(icon)

        btn.setIconSize(QSize(32, 32))
        btn.setMinimumSize(QSize(40, 40))

        btn.setToolTip(alt_text)

        self.bar_layout.addWidget(btn)

        return btn

    def _add_toggle_button(self, alt_text: str, icon: str | None = None, value: bool = True):
        """ Add a toggle button and wire it up"""
        btn = self._toolbar_button(alt_text, icon)

        btn.setCheckable(True)
        btn.setChecked(value)

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

        self.prettify = self._add_toggle_button("Use cosmetic representation instead of direct input",
                                                icon="cosmetic",
                                                value=settings.prettify)


        self.bar_layout.addSpacerItem(QSpacerItem(20, 0, QSizePolicy.Minimum, QSizePolicy.Minimum))

        #
        # Various visual properties
        #

        self.atom_or_moments = self._add_toggle_button("Switch between showing atoms and moments",
                                                       icon="momentatoms",
                                                       value=settings.show_atoms_not_moments)




        self.moment_scale_slider = self._add_slider(
            0, 1, settings.atom_moment_scaling,
            left_label=IconWidget("small_moments", "Smaller Sites"),
            right_label=IconWidget("big_moments", "Larger Sites"),
            alt_text="Site scale factor")


        self.show_nonmagnetic = self._add_toggle_button("Show non-magnetic sites",
                                                        icon="nonmagnetic",
                                                        value=settings.show_nonmagnetic_atoms)

        self.scale_atoms = self._add_toggle_button("Scale atoms by atomic radii",
                                                   icon="atomsizes",
                                                   value=settings.use_atomic_radii)

        self.coupling_scale_slider = self._add_slider(
            0, 1, settings.coupling_scaling,
            left_label=IconWidget("small_exchange", "Smaller Exchanges"),
            right_label=IconWidget("big_exchange", "Larger Exchanges"),
            alt_text="Exchange thickness")



        self.bar_layout.addSpacerItem(QSpacerItem(20, 0, QSizePolicy.Minimum, QSizePolicy.Minimum))

        self.show_cartesian_axes = self._add_toggle_button("Show Cartesian axes",
                                                           icon="cartesianaxes",
                                                           value=settings.show_cartesian_axes)

        self.show_lattice_axes = self._add_toggle_button("Show lattice aligned axes",
                                                         icon="latticeaxes",
                                                         value=settings.show_lattice_axes)

        self.orthogonal_lattice_axes = self._add_toggle_button("Show orthogonal lattice axes",
                                                               icon="latticeaxesorth",
                                                               value=settings.orthogonal_lattice_axes)

        self.reset_view = self._toolbar_button("Reset view", icon="datum")

        self.reset_view.clicked.connect(self.on_reset_view_clicked)

        # Pad right, and set layout
        self.bar_layout.addSpacerItem(QSpacerItem(0, 0, QSizePolicy.Expanding, QSizePolicy.Minimum))
        self.setLayout(self.bar_layout)

        # Scaling properties
        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Minimum)
        self.bar_layout.setContentsMargins(0, 0, 0, 0)
        self.bar_layout.setSpacing(0)

        # TODO: Temporally disabled until implemented
        self.show_anisotropies.setVisible(False)
        self.atom_or_moments.setVisible(False)
        self.scale_atoms.setVisible(False)

    def _on_change(self):
        """ Called when anything changes, send signal """
        self.displayOptionsChanged.emit()

    def on_reset_view_clicked(self):
        """ Called when the view reset button is pressed, send a signal"""
        self.requestViewReset.emit()

    def display_options(self) -> DisplayOptions:
        """ Get the current render options """
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
            show_cartesian_axes = self.show_cartesian_axes.isChecked(),
            show_lattice_axes = self.show_lattice_axes.isChecked(),
            orthogonal_lattice_axes = self.orthogonal_lattice_axes.isChecked(),
            prettify = self.prettify.isChecked())

    def save_settings(self):
        """ Save the current settings"""
        with open(self.settings_filename, 'w') as file:
            file.write(self.display_options().serialise())

    def closeEvent(self, event):
        """ Qt override, window closed"""
        # we want to save the settings before exiting
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
