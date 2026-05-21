""" Window that shows when settings button is clicked """

from dataclasses import replace

from PySide6.QtWidgets import QDialog, QVBoxLayout, QLabel, QDialogButtonBox, QColorDialog, QWidget, QGridLayout

from pyspinw.gui.displayoptions import DisplayOptions
from pyspinw.gui.util import QRightLabel, QColorPatch


class SettingsDialog(QDialog):
    """ Dialog with settings options """

    def __init__(self, current_options: DisplayOptions, parent=None):
        self._current_options = current_options
        super().__init__(parent=parent)

        self.setWindowTitle("Settings")

        layout = QVBoxLayout(self)

        inner_widget = QWidget(self)
        inner_layout = QGridLayout()
        inner_widget.setLayout(inner_layout)

        self._background = QColorPatch(*current_options.background_color, parent=self)
        inner_layout.addWidget(QRightLabel("Background Color"), 0, 0)
        inner_layout.addWidget(self._background, 0, 1)

        self._sites = QColorPatch(*current_options.default_site_color, parent=self)
        inner_layout.addWidget(QRightLabel("Sites Default Color"), 1, 0)
        inner_layout.addWidget(self._sites, 1, 1)

        self._exchanges = QColorPatch(*current_options.default_exchange_color, parent=self)
        inner_layout.addWidget(QRightLabel("Exchange Default Color"), 2, 0)
        inner_layout.addWidget(self._exchanges, 2, 1)

        # Add inner widget
        layout.addWidget(inner_widget)

        # OK / Cancel buttons
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok |
            QDialogButtonBox.StandardButton.Cancel
        )

        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout.addWidget(buttons)

    def new_display_options(self):
        """ Get the updated display options """
        return replace(
            self._current_options,
            background_color = self._background.rgb,
            default_exchange_color = self._exchanges.rgb,
            default_site_color = self._sites.rgb)


if __name__ == "__main__":
    options = DisplayOptions()

    from PySide6.QtWidgets import QApplication

    app = QApplication()
    widget = SettingsDialog(options)

    if widget.exec() == QDialog.Accepted:
        value = widget.new_display_options()
        print(value)
