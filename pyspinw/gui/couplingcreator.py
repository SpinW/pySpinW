from dataclasses import dataclass

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QWidget, QHBoxLayout, QLineEdit, QVBoxLayout, QPushButton, QApplication, QComboBox, \
    QSpacerItem, QSizePolicy, QLabel

from pyspinw.batch_couplings import batch_couplings, default_naming_pattern
from pyspinw.calculations.spinwave import Coupling
from pyspinw.coupling import couplings as coupling_classes
from pyspinw.coupling import coupling_lookup
from pyspinw.gui.couplingtable import CouplingTable
from pyspinw.gui.helperwidgets.couplingtypecombo import CouplingTypeCombo
from pyspinw.gui.helperwidgets.floatfield import FloatField
from pyspinw.gui.helperwidgets.parametertable import ParameterTable, ParameterExists
from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell

@dataclass
class CreationParameters:
    format_string: str
    coupling_type: type[Coupling]
    min_distance: float
    max_distance: float
    max_order: int | None

class CouplingCreator(QWidget):

    ok_clicked = Signal()
    cancel_clicked = Signal()

    def __init__(self, sites: list[LatticeSite], unit_cell: UnitCell, parent=None):
        super().__init__(parent=parent)

        self.sites = sites
        self.unit_cell = unit_cell
        self.couplings = None

        # Table for viewing output
        self.output_table = CouplingTable(editable=False)

        # Label for couplings
        self.coupling_label = QLabel("Couplings (0 found):")

        # Overall layout

        field_widget = ParameterTable()

        button_layout = QHBoxLayout()
        button_widget = QWidget()
        button_widget.setLayout(button_layout)

        main_layout = QVBoxLayout()
        main_layout.addWidget(field_widget)
        main_layout.addWidget(self.coupling_label)
        main_layout.addWidget(self.output_table, stretch=1)
        main_layout.addWidget(button_widget)

        self.setLayout(main_layout)

        # Fields

        name_string = QLineEdit(default_naming_pattern)
        field_widget.add_parameter(name="name_string",
                                   widget=name_string,
                                   value_getter=name_string.text,
                                   change_signal=name_string.textChanged,
                                   display_alias="Formatting")

        min_distance = FloatField(0, bottom=0, slider_top=3*unit_cell.main_diagonal_length)
        field_widget.add_parameter(name="min_distance",
                                   widget=min_distance,
                                   value_getter=lambda: min_distance.value,
                                   change_signal=min_distance.changed,
                                   display_alias="Min. Distance")

        max_distance = FloatField(int(50*unit_cell.main_diagonal_length)/100, bottom=0, slider_top=3*unit_cell.main_diagonal_length)
        field_widget.add_parameter(name="max_distance",
                                   widget=max_distance,
                                   value_getter=lambda: max_distance.value,
                                   change_signal=max_distance.changed,
                                   display_alias="Max. Distance")

        # Maximum order combo
        max_order = QComboBox()
        max_order.addItem("None")
        for i in range(1, 6):
            max_order.addItem(str(i))
        max_order.setCurrentText("None")
        field_widget.add_parameter(name="max_order",
                                   widget=max_order,
                                   value_getter=max_order.currentText,
                                   change_signal=max_order.currentTextChanged,
                                   display_alias="Max. Order")

        # Coupling type combo
        coupling_type = CouplingTypeCombo(coupling_classes[0].coupling_type)
        field_widget.add_parameter(name="coupling_type",
                                   widget=coupling_type,
                                   value_getter=coupling_type.currentText,
                                   change_signal=coupling_type.currentTextChanged,
                                   display_alias="Coupling Type")

        coupling_type.currentTextChanged.connect(self._update_fields)

        # Add fields for all the different kinds of couplings
        for coupling in coupling_classes:
            for parameter, default in zip(coupling.parameters, coupling.parameter_defaults):
                field = FloatField(default, slider_bottom=-50, slider_top=50)

                try:
                    field_widget.add_parameter(name=parameter,
                                               widget=field,
                                               value_getter=field.get_value,
                                               change_signal=field.changed)


                except ParameterExists: # Very specific error for duplicated keys, we should ignore it in this case
                    pass

        field_widget.set_all_visible(False)
        self._always_shown = ["name_string", "min_distance", "max_distance", "max_order", "coupling_type"]
        self.field_widget = field_widget

        self.field_widget.changed.connect(self._on_field_changed)

        # Buttons

        self.cancel_button = QPushButton("Cancel")
        self.ok_button = QPushButton("OK")

        button_layout.addWidget(self.cancel_button)
        button_layout.addSpacerItem(QSpacerItem(20,20, QSizePolicy.Expanding, QSizePolicy.Minimum))
        button_layout.addWidget(self.ok_button)

        self.cancel_button.clicked.connect(self._on_cancel_clicked)
        self.ok_button.clicked.connect(self._on_ok_clicked)

        # Set state

        self._update_fields()
        self._update()

    def _current_type(self):
        return coupling_lookup[self.field_widget.get_value("coupling_type")]

    def _parameters(self) -> CreationParameters:
        coupling_type = self._current_type()
        max_order_string = self.field_widget.get_value("max_order")
        max_order = None if max_order_string == "None" else int(max_order_string)

        return CreationParameters(
            format_string=self.field_widget.get_value("name_string"),
            coupling_type=coupling_type,
            min_distance=self.field_widget.get_value("min_distance"),
            max_distance=self.field_widget.get_value("max_distance"),
            max_order=max_order)

    def _update_label(self):
        self.coupling_label.setText("Couplings (%i found)" % len(self.couplings))

    def _calculate_couplings(self) -> list[Coupling]:
        parameters = self._parameters()
        coupling_type = self._current_type()

        abstract_couplings = batch_couplings(
            sites=self.sites,
            unit_cell=self.unit_cell,
            max_distance=parameters.max_distance,
            naming_pattern=parameters.format_string,
            type_symbol=coupling_type.short_string)

        parameter_dict = {parameter: self.field_widget.get_value(parameter)
                          for parameter in coupling_type.parameters}

        kept_couplings = []
        for coupling in abstract_couplings:
            if parameters.max_order is not None and coupling.order > parameters.max_order:
                continue

            if parameters.min_distance <= coupling.distance <= parameters.max_distance:

                kept_couplings.append(coupling_type(
                    name=coupling.name,
                    site_1=coupling.site_1,
                    site_2=coupling.site_2,
                    cell_offset=coupling.cell_offset,
                    **parameter_dict))

        return kept_couplings

    def _update_fields(self):
        self.field_widget.set_all_visible(False)

        current_type = self._current_type()

        for parameter in self._always_shown:
            self.field_widget.set_parameter_visible(parameter, True)

        for parameter in current_type.parameters:
            self.field_widget.set_parameter_visible(parameter, True)

    def _on_field_changed(self):
        self._update()

    def _update(self):
        self.couplings = self._calculate_couplings()
        self.output_table.couplings = self.couplings
        self._update_label()

    def _on_ok_clicked(self):
        self.ok_clicked.emit()

    def _on_cancel_clicked(self):
        self.cancel_clicked.emit()



if __name__ == "__main__":
    app = QApplication([])

    sites = [LatticeSite.create(0,0,0,0,0,1, "A"),
             LatticeSite.create(0.1,0,0,0,0,1, "B")]

    unit_cell = UnitCell(1,1,1)

    coupling_creator = CouplingCreator(sites, unit_cell)
    coupling_creator.show()

    app.exec_()
