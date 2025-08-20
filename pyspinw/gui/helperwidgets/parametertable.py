""" A table for editing different kinds of parameters, it can change dynamically"""

from collections import defaultdict
from typing import Callable, Any

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QWidget, QGridLayout

from pyspinw.gui.helperwidgets.misc import QRightLabel, QLeftLabel

class ParameterExists(Exception):
    """ Exception thrown when trying to display two parameters with the same name"""

class ParameterTable(QWidget):
    """ A dynamic table for showing different kinds of parameters in a nice grid"""

    changed = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)

        self._row_count = 0
        self._name_row_lookup = {}

        self._name_value_functions: dict[str, Callable[[], Any]] = {}
        self._minimum_heights = defaultdict[str](int)

        self.grid_layout = QGridLayout()
        self.setLayout(self.grid_layout)

    def set_all_visible(self, visible: bool):
        """ Show/hide everything"""
        for row in range(self._row_count):
            self._set_row_visible(row, visible)

    def set_parameter_visible(self, name: str, visible: bool):
        """ Show/hide a specific name """
        self._set_row_visible(self._name_row_lookup[name], visible)


    def add_parameter(self,
                      name: str,
                      widget: QWidget,
                      value_getter: Callable[[], Any],
                      change_signal: Signal,
                      units: str | None = None,
                      display_alias: str | None = None):
        """Add a parameter to the grid

        :param name: Name to use for this parameter
        :param widget: Widget to put in the table
        :param value_getter: function to get the value from the specified parameter, typing is difficult here,
                             so it is left to the user to fix
        :param change_signal: signal that the added widget should emit from on change
        :param units: units for unit field
        :param display_alias: Name to show instead of the parameter name
        """
        if name in self._name_row_lookup:
            raise ParameterExists(f"Parameter {name} already exists")

        display_name = name if display_alias is None else display_alias

        self.grid_layout.addWidget(QRightLabel(display_name), self._row_count, 0)
        self.grid_layout.addWidget(widget, self._row_count, 1)
        self._name_row_lookup[name] = self._row_count
        self._name_value_functions[name] = value_getter

        change_signal.connect(self._on_change)

        if units is not None:
            self.grid_layout.addWidget(QLeftLabel(units), self._row_count, 2)

        self._row_count += 1

    def _set_row_visible(self, row_index, visible_state):
        """ Show or hide a row, pretty messy (hence having a dedicated class to handle all this)"""
        if visible_state:

            # Use current height if it exists
            current_height = self.grid_layout.rowMinimumHeight(row_index)
            if current_height == 0:
                target_height = self._minimum_heights[row_index]
            else:
                target_height = current_height

            # Show all the widgets in the row
            for col in range(self.grid_layout.columnCount()):
                item = self.grid_layout.itemAtPosition(row_index, col)
                if item:
                    widget = item.widget()
                    if widget:
                        widget.show()

            # Reembiggen row
            self.grid_layout.setRowMinimumHeight(row_index, target_height)
            self.grid_layout.setRowStretch(row_index, 1)

        else:
            # Save the height
            self._minimum_heights[row_index] = self.grid_layout.rowMinimumHeight(row_index)

            # Hide all the widgets
            for col in range(self.grid_layout.columnCount()):
                item = self.grid_layout.itemAtPosition(row_index, col)
                if item:
                    widget = item.widget()
                    if widget:
                        widget.hide()

            # Disembiggen row
            self.grid_layout.setRowMinimumHeight(row_index, 0)
            self.grid_layout.setRowStretch(row_index, 0)


    def get_value(self, parameter_name: str):
        """ Get the value of a given parameter

        Note that this is not typed as different parameters can have different kinds of values
        """
        return self._name_value_functions[parameter_name]()

    def _on_change(self):
        self.changed.emit()
