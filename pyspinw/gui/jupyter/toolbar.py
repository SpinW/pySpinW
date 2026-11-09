"""Jupyter ipywidgets toolbar and settings panel for the crystal viewer."""

import base64
from collections.abc import Callable
from dataclasses import replace
from importlib import resources
import os

import ipywidgets as widgets

from pyspinw.gui.displayoptions import DisplayOptions


def rgb_to_hex(rgb: tuple[float, float, float]) -> str:
    """Convert an RGB float tuple in ``[0, 1]`` to a ``#rrggbb`` hex string.

    Parameters
    ----------
    rgb : tuple of (float, float, float)
        RGB values in the range ``[0, 1]``.

    Returns
    -------
    str
        Hex color string in ``#rrggbb`` format.
    """
    r = int(max(0, min(255, round(rgb[0] * 255))))
    g = int(max(0, min(255, round(rgb[1] * 255))))
    b = int(max(0, min(255, round(rgb[2] * 255))))
    return f"#{r:02x}{g:02x}{b:02x}"


def hex_to_rgb(hex_color: str) -> tuple[float, float, float]:
    """Convert a ``#rrggbb`` or ``#rgb`` hex color string to an RGB float tuple in ``[0, 1]``.

    Parameters
    ----------
    hex_color : str
        Hex color string starting with ``#``.

    Returns
    -------
    tuple of (float, float, float)
        RGB values in the range ``[0, 1]``.
    """
    value = hex_color.lstrip("#")
    if len(value) == 3:
        value = "".join(ch * 2 for ch in value)
    r = int(value[0:2], 16) / 255.0
    g = int(value[2:4], 16) / 255.0
    b = int(value[4:6], 16) / 255.0
    return (r, g, b)


def load_svg_data_uri(icon_name: str, dark: bool = False) -> str:
    """Load an SVG icon from ``pyspinw.gui.icons.svg`` and encode it as a ``data:`` URI.

    Parameters
    ----------
    icon_name : str
        Base name of the SVG icon without extension.
    dark : bool, optional
        Whether to use the ``-dark`` icon variant. Default is ``False``.

    Returns
    -------
    str
        Base64-encoded ``data:image/svg+xml;base64,...`` URI string.
    """
    variant = f"{icon_name}-dark.svg" if dark else f"{icon_name}.svg"
    svg_files = resources.files("pyspinw.gui.icons.svg")
    target = svg_files / variant
    if not target.is_file():
        target = svg_files / f"{icon_name}.svg"
    raw_bytes = target.read_bytes()
    b64 = base64.b64encode(raw_bytes).decode("ascii")
    return f"data:image/svg+xml;base64,{b64}"


def svg_img_html(icon_name: str, alt_text: str = "", size: int = 22, dark: bool = False) -> str:
    """Build an HTML ``<img>`` tag displaying an SVG icon from ``pyspinw.gui.icons.svg``.

    Parameters
    ----------
    icon_name : str
        Base name of the SVG icon without extension.
    alt_text : str, optional
        Tooltip and alt text for the image. Default is ``""``.
    size : int, optional
        Width and height in pixels. Default is ``22``.
    dark : bool, optional
        Whether to use the ``-dark`` icon variant. Default is ``False``.

    Returns
    -------
    str
        HTML string containing the ``<img>`` element.
    """
    uri = load_svg_data_uri(icon_name, dark=dark)
    escaped_alt = alt_text.replace('"', "&quot;")
    return (
        f'<img src="{uri}" width="{size}" height="{size}" '
        f'title="{escaped_alt}" alt="{escaped_alt}" '
        f'style="vertical-align: middle; display: inline-block;" />'
    )


class JupyterIconWidget(widgets.HTML):
    """Display an SVG icon inside an ``ipywidgets.HTML`` widget.

    Parameters
    ----------
    icon_name : str
        Name of the SVG icon to load.
    alt_text : str, optional
        Tooltip text for the widget. Default is ``""``.
    size : int, optional
        Icon size in pixels. Default is ``22``.
    """

    def __init__(self, icon_name: str, alt_text: str = "", size: int = 22):
        super().__init__(
            value=svg_img_html(icon_name, alt_text=alt_text, size=size),
            layout=widgets.Layout(
                width=f"{size + 4}px",
                height=f"{size + 4}px",
                margin="0px 2px",
                display="flex",
                align_items="center",
                justify_content="center",
            ),
        )
        self.icon_name = icon_name
        self.alt_text = alt_text


class JupyterGraphicalSlider(widgets.HBox):
    """Slider widget with left and right icon labels for the Jupyter toolbar.

    Parameters
    ----------
    min_value : float
        Minimum value represented by the slider.
    max_value : float
        Maximum value represented by the slider.
    start_value : float
        Initial value of the slider.
    left_label : widgets.Widget
        Widget displayed to the left of the slider.
    right_label : widgets.Widget
        Widget displayed to the right of the slider.
    alt_text : str
        Tooltip text for the slider and its labels.
    """

    def __init__(
        self,
        min_value: float,
        max_value: float,
        start_value: float,
        left_label: widgets.Widget,
        right_label: widgets.Widget,
        alt_text: str,
    ):
        self.min_value = float(min_value)
        self.max_value = float(max_value)
        clamped_start = max(self.min_value, min(self.max_value, float(start_value)))
        step = (self.max_value - self.min_value) / 1000.0 if self.max_value > self.min_value else 0.001

        self.slider = widgets.FloatSlider(
            value=clamped_start,
            min=self.min_value,
            max=self.max_value,
            step=step,
            readout=False,
            tooltip=alt_text,
            layout=widgets.Layout(width="110px", margin="0px 2px"),
        )

        super().__init__(
            children=[left_label, self.slider, right_label],
            layout=widgets.Layout(align_items="center", margin="0px 6px"),
        )

        self._callbacks: list[Callable[[], None]] = []
        self.slider.observe(self._on_slider_change, names="value")

    def _on_slider_change(self, _change):
        for cb in self._callbacks:
            cb()

    def on_value_changed(self, callback: Callable[[], None]):
        """Register a callback invoked when the slider value changes.

        Parameters
        ----------
        callback : callable
            Zero-argument function called when the slider moves.
        """
        self._callbacks.append(callback)

    def value(self) -> float:
        """Return the current float value of the slider.

        Returns
        -------
        float
            Current slider value in ``[min_value, max_value]``.
        """
        return float(self.slider.value)


class JupyterSettingsPanel(widgets.VBox):
    """Collapsible color settings panel replicating ``SettingsDialog`` in ``ipywidgets``.

    Parameters
    ----------
    current_options : DisplayOptions
        Initial display options used to populate the color pickers.
    """

    def __init__(self, current_options: DisplayOptions):
        self._current_options = current_options
        self._apply_callbacks: list[Callable[[DisplayOptions], None]] = []
        self._cancel_callbacks: list[Callable[[], None]] = []

        picker_layout = widgets.Layout(width="260px")

        self.background_picker = widgets.ColorPicker(
            description="Background:",
            value=rgb_to_hex(current_options.background_color),
            concise=False,
            style={"description_width": "130px"},
            layout=picker_layout,
        )
        self.sites_picker = widgets.ColorPicker(
            description="Sites Default:",
            value=rgb_to_hex(current_options.default_site_color),
            concise=False,
            style={"description_width": "130px"},
            layout=picker_layout,
        )
        self.exchanges_picker = widgets.ColorPicker(
            description="Exchange Default:",
            value=rgb_to_hex(current_options.default_exchange_color),
            concise=False,
            style={"description_width": "130px"},
            layout=picker_layout,
        )
        self.selected_picker = widgets.ColorPicker(
            description="Selection:",
            value=rgb_to_hex(current_options.selected_color),
            concise=False,
            style={"description_width": "130px"},
            layout=picker_layout,
        )
        self.hover_picker = widgets.ColorPicker(
            description="Hover:",
            value=rgb_to_hex(current_options.hover_color),
            concise=False,
            style={"description_width": "130px"},
            layout=picker_layout,
        )
        self.selected_hover_picker = widgets.ColorPicker(
            description="Selected Hover:",
            value=rgb_to_hex(current_options.selected_hover_color),
            concise=False,
            style={"description_width": "130px"},
            layout=picker_layout,
        )
        self.anisotropies_picker = widgets.ColorPicker(
            description="Anisotropies:",
            value=rgb_to_hex(current_options.anisotropy_color),
            concise=False,
            style={"description_width": "130px"},
            layout=picker_layout,
        )

        self.ok_button = widgets.Button(
            description="OK",
            button_style="primary",
            layout=widgets.Layout(width="80px"),
        )
        self.cancel_button = widgets.Button(
            description="Cancel",
            layout=widgets.Layout(width="80px"),
        )

        self.ok_button.on_click(self._on_ok_clicked)
        self.cancel_button.on_click(self._on_cancel_clicked)

        grid = widgets.HBox(
            [
                widgets.VBox(
                    [
                        self.background_picker,
                        self.sites_picker,
                        self.exchanges_picker,
                        self.anisotropies_picker,
                    ]
                ),
                widgets.VBox(
                    [
                        self.selected_picker,
                        self.hover_picker,
                        self.selected_hover_picker,
                    ]
                ),
            ]
        )
        buttons_row = widgets.HBox([self.ok_button, self.cancel_button], layout=widgets.Layout(margin="4px 0px"))

        super().__init__(
            children=[grid, buttons_row],
            layout=widgets.Layout(
                display="none",
                border="1px solid #888",
                padding="6px",
                margin="4px 0px",
            ),
        )

    def sync_from_options(self, options: DisplayOptions):
        """Update all color pickers to match ``options``.

        Parameters
        ----------
        options : DisplayOptions
            Display options to read colors from.
        """
        self._current_options = options
        self.background_picker.value = rgb_to_hex(options.background_color)
        self.sites_picker.value = rgb_to_hex(options.default_site_color)
        self.exchanges_picker.value = rgb_to_hex(options.default_exchange_color)
        self.selected_picker.value = rgb_to_hex(options.selected_color)
        self.hover_picker.value = rgb_to_hex(options.hover_color)
        self.selected_hover_picker.value = rgb_to_hex(options.selected_hover_color)
        self.anisotropies_picker.value = rgb_to_hex(options.anisotropy_color)

    def new_display_options(self, base_options: DisplayOptions | None = None) -> DisplayOptions:
        """Return a ``DisplayOptions`` instance incorporating the current picker colors.

        Parameters
        ----------
        base_options : DisplayOptions or None, optional
            Options to replace color fields on. Defaults to ``self._current_options``.

        Returns
        -------
        DisplayOptions
            Updated display options.
        """
        target = self._current_options if base_options is None else base_options
        return replace(
            target,
            background_color=hex_to_rgb(self.background_picker.value),
            default_site_color=hex_to_rgb(self.sites_picker.value),
            default_exchange_color=hex_to_rgb(self.exchanges_picker.value),
            selected_color=hex_to_rgb(self.selected_picker.value),
            hover_color=hex_to_rgb(self.hover_picker.value),
            selected_hover_color=hex_to_rgb(self.selected_hover_picker.value),
            anisotropy_color=hex_to_rgb(self.anisotropies_picker.value),
        )

    def on_apply(self, callback: Callable[[DisplayOptions], None]):
        """Register a callback invoked when the OK button is clicked.

        Parameters
        ----------
        callback : callable
            Function accepting the updated ``DisplayOptions``.
        """
        self._apply_callbacks.append(callback)

    def on_cancel(self, callback: Callable[[], None]):
        """Register a callback invoked when the Cancel button is clicked.

        Parameters
        ----------
        callback : callable
            Zero-argument function.
        """
        self._cancel_callbacks.append(callback)

    def _on_ok_clicked(self, _btn):
        updated = self.new_display_options()
        self.layout.display = "none"
        for cb in self._apply_callbacks:
            cb(updated)

    def _on_cancel_clicked(self, _btn):
        self.sync_from_options(self._current_options)
        self.layout.display = "none"
        for cb in self._cancel_callbacks:
            cb()


class JupyterToolbar(widgets.VBox):
    """Jupyter ``ipywidgets`` toolbar replicating ``DisplayOptionsToolbar``.

    Parameters
    ----------
    initial_display_options : DisplayOptions or None, optional
        Initial display options to populate the toolbar. If ``None``, loads
        from ``render_settings.conf`` or uses defaults. Default is ``None``.
    max_size_object_scaling : float, optional
        Maximum value for the object scaling sliders. Default is ``10.0``.
    save_config_on_exit : bool, optional
        Whether to save display settings when updated. Default is ``True``.

    Raises
    ------
    TypeError
        If ``initial_display_options`` is provided and is not an instance of
        ``DisplayOptions``.
    """

    settings_filename = "render_settings.conf"

    def __init__(
        self,
        initial_display_options: DisplayOptions | None = None,
        max_size_object_scaling: float = 10.0,
        save_config_on_exit: bool = True,
    ):
        self.save_config_on_exit = save_config_on_exit

        if initial_display_options is None:
            if os.path.exists(self.settings_filename):
                with open(self.settings_filename, "r", encoding="utf-8") as file:
                    settings = DisplayOptions.deserialise(file.read())
            else:
                settings = DisplayOptions()
        else:
            if not isinstance(initial_display_options, DisplayOptions):
                raise TypeError("Expected initial_display_options to be an instance of DisplayOptions")
            settings = initial_display_options

        self._background_color = settings.background_color
        self._sites_color = settings.default_site_color
        self._exchanges_color = settings.default_exchange_color
        self._selected_color = settings.selected_color
        self._hover_color = settings.hover_color
        self._selected_hover_color = settings.selected_hover_color
        self._anisotropy_color = settings.anisotropy_color

        self._display_options_callbacks: list[Callable[[], None]] = []
        self._reset_view_callbacks: list[Callable[[], None]] = []
        self._snapshot_callbacks: list[Callable[[], None]] = []

        self.bar_1_items: list[widgets.Widget] = []
        self.bar_2_items: list[widgets.Widget] = []

        # Top Bar Controls
        self.prettify = self._add_toggle_button(
            self.bar_1_items,
            alt_text="Use cosmetic representation instead of direct input",
            icon="cosmetic",
            value=settings.prettify,
        )

        self.atom_or_spins = self._add_toggle_button(
            self.bar_1_items,
            alt_text="Switch between showing atoms and spins",
            icon="momentatoms",
            value=settings.show_spins_not_atoms,
        )
        self.show_nonmagnetic = self._add_toggle_button(
            self.bar_1_items,
            alt_text="Show non-magnetic sites",
            icon="nonmagnetic",
            value=settings.show_nonmagnetic_atoms,
        )
        self.scale_atoms = self._add_toggle_button(
            self.bar_1_items,
            alt_text="Scale atoms by atomic radii",
            icon="atomsizes",
            value=settings.use_atomic_radii,
            visible=False,
        )

        self.show_cartesian_axes = self._add_toggle_button(
            self.bar_1_items,
            alt_text="Show Cartesian axes",
            icon="cartesianaxes",
            value=settings.show_cartesian_axes,
        )
        self.show_lattice_axes = self._add_toggle_button(
            self.bar_1_items,
            alt_text="Show lattice aligned axes",
            icon="latticeaxes",
            value=settings.show_lattice_axes,
        )
        self.orthogonal_lattice_axes = self._add_toggle_button(
            self.bar_1_items,
            alt_text="Show orthogonal lattice axes",
            icon="latticeaxesorth",
            value=settings.orthogonal_lattice_axes,
        )
        self.orthographic = self._add_toggle_button(
            self.bar_1_items,
            alt_text="Toggle Orthographic/Perspective",
            icon="perspective",
            value=settings.perspective,
        )

        self.reset_view = self._toolbar_button(
            self.bar_1_items,
            alt_text="Reset view",
            icon="datum",
        )
        self.reset_view.on_click(self._on_reset_view_clicked)

        self.snapshot_button = self._toolbar_button(
            self.bar_1_items,
            alt_text="Snapshot",
            icon="camera",
        )
        self.snapshot_button.on_click(self._on_snapshot_clicked)

        self.settings_button = self._toolbar_button(
            self.bar_1_items,
            alt_text="Settings",
            icon="settings",
        )
        self.settings_button.on_click(self._on_settings_clicked)

        # Lower Bar Controls
        self.show_cell = self._add_toggle_button(
            self.bar_2_items,
            alt_text="Show Unit Cell",
            icon="unitcell",
            value=settings.show_unit_cell,
        )
        self.show_supercell = self._add_toggle_button(
            self.bar_2_items,
            alt_text="Show Supercell",
            icon="supercell",
            value=settings.show_supercell,
        )
        self.show_sites = self._add_toggle_button(
            self.bar_2_items,
            alt_text="Show Sites",
            icon="show_moments",
            value=settings.show_sites,
        )
        self.show_exchanges = self._add_toggle_button(
            self.bar_2_items,
            alt_text="Show Exchanges",
            icon="show_exchanges",
            value=settings.show_exchanges,
        )
        self.show_anisotropies = self._add_toggle_button(
            self.bar_2_items,
            alt_text="Show Anisotropies",
            icon="anisotropies",
            value=settings.show_anisotropies,
        )

        self.spin_scale_slider = self._add_slider(
            self.bar_2_items,
            0.0,
            max_size_object_scaling,
            settings.atom_spin_scaling,
            left_label=JupyterIconWidget("small_moments", "Smaller Sites"),
            right_label=JupyterIconWidget("big_moments", "Larger Sites"),
            alt_text="Site scale factor",
        )
        self.exchange_scale_slider = self._add_slider(
            self.bar_2_items,
            0.0,
            max_size_object_scaling,
            settings.exchange_scaling,
            left_label=JupyterIconWidget("small_exchange", "Smaller Exchanges"),
            right_label=JupyterIconWidget("big_exchange", "Larger Exchanges"),
            alt_text="Exchange thickness",
        )
        self.anisotropy_scale_slider = self._add_slider(
            self.bar_2_items,
            0.0,
            max_size_object_scaling,
            settings.anisotropy_scaling,
            left_label=JupyterIconWidget("small_anisotropies", "Smaller Bubbles"),
            right_label=JupyterIconWidget("big_anisotropies", "Larger Bubbles"),
            alt_text="Anisotropy visualisation size",
        )

        self.bar_1 = widgets.HBox(
            self.bar_1_items,
            layout=widgets.Layout(flex_flow="row wrap", align_items="center", margin="2px 0px"),
        )
        self.bar_2 = widgets.HBox(
            self.bar_2_items,
            layout=widgets.Layout(flex_flow="row wrap", align_items="center", margin="2px 0px"),
        )

        self.settings_panel = JupyterSettingsPanel(self.display_options())
        self.settings_panel.on_apply(self._on_settings_applied)

        super().__init__(
            children=[self.bar_1, self.bar_2, self.settings_panel],
            layout=widgets.Layout(width="100%", margin="0px 0px 4px 0px"),
        )

    def _toolbar_button(
        self,
        target_list: list[widgets.Widget],
        alt_text: str,
        icon: str | None = None,
    ) -> widgets.Button:
        """Create an action button paired with an SVG icon tile.

        Parameters
        ----------
        target_list : list of widgets.Widget
            List of toolbar widgets to append the button tile to.
        alt_text : str
            Tooltip text for the button.
        icon : str or None, optional
            SVG icon name from ``pyspinw.gui.icons.svg``. Default is ``None``.

        Returns
        -------
        widgets.Button
            Created button widget.
        """
        btn = widgets.Button(
            description="",
            tooltip=alt_text,
            layout=widgets.Layout(width="36px", height="32px", margin="0px"),
        )
        if icon is not None:
            icon_widget = JupyterIconWidget(icon, alt_text=alt_text, size=20)
            tile = widgets.HBox(
                [icon_widget, btn],
                layout=widgets.Layout(
                    align_items="center",
                    border="1px solid #555",
                    margin="1px 2px",
                    padding="1px 3px",
                ),
            )
            target_list.append(tile)
        else:
            target_list.append(btn)
        return btn

    def _add_toggle_button(
        self,
        target_list: list[widgets.Widget],
        alt_text: str,
        icon: str | None = None,
        value: bool = True,
        visible: bool = True,
    ) -> widgets.ToggleButton:
        """Create a checkable ``ToggleButton`` paired with an SVG icon tile.

        Parameters
        ----------
        target_list : list of widgets.Widget
            List of toolbar widgets to append the toggle tile to.
        alt_text : str
            Tooltip text for the toggle button.
        icon : str or None, optional
            SVG icon name from ``pyspinw.gui.icons.svg``. Default is ``None``.
        value : bool, optional
            Initial checked state. Default is ``True``.
        visible : bool, optional
            Whether the control is visible in the toolbar. Default is ``True``.

        Returns
        -------
        widgets.ToggleButton
            Created toggle button widget.
        """
        btn = widgets.ToggleButton(
            value=bool(value),
            description="",
            tooltip=alt_text,
            button_style="info" if value else "",
            layout=widgets.Layout(width="28px", height="28px", margin="0px"),
        )

        def _sync_style(change):
            btn.button_style = "info" if change["new"] else ""
            self._on_change()

        btn.observe(_sync_style, names="value")

        if icon is not None:
            icon_widget = JupyterIconWidget(icon, alt_text=alt_text, size=20)
            tile = widgets.HBox(
                [icon_widget, btn],
                layout=widgets.Layout(
                    display="flex" if visible else "none",
                    align_items="center",
                    border="1px solid #555",
                    margin="1px 2px",
                    padding="1px 3px",
                ),
            )
            target_list.append(tile)
        else:
            if not visible:
                btn.layout.display = "none"
            target_list.append(btn)

        return btn

    def _add_slider(
        self,
        target_list: list[widgets.Widget],
        min_value: float,
        max_value: float,
        start_value: float,
        left_label: widgets.Widget,
        right_label: widgets.Widget,
        alt_text: str,
    ) -> JupyterGraphicalSlider:
        """Create a ``JupyterGraphicalSlider`` and wire it to toolbar change notifications.

        Parameters
        ----------
        target_list : list of widgets.Widget
            List of toolbar widgets to append the slider to.
        min_value : float
            Minimum slider value.
        max_value : float
            Maximum slider value.
        start_value : float
            Initial slider value.
        left_label : widgets.Widget
            Widget displayed on the left of the slider.
        right_label : widgets.Widget
            Widget displayed on the right of the slider.
        alt_text : str
            Tooltip text for the slider.

        Returns
        -------
        JupyterGraphicalSlider
            Created graphical slider widget.
        """
        slider = JupyterGraphicalSlider(
            min_value=min_value,
            max_value=max_value,
            start_value=start_value,
            left_label=left_label,
            right_label=right_label,
            alt_text=alt_text,
        )
        slider.on_value_changed(self._on_change)
        target_list.append(slider)
        return slider

    def on_display_options_changed(self, callback: Callable[[], None]):
        """Register a callback invoked whenever any display option control changes.

        Parameters
        ----------
        callback : callable
            Zero-argument function.
        """
        self._display_options_callbacks.append(callback)

    def on_request_view_reset(self, callback: Callable[[], None]):
        """Register a callback invoked when the Reset View button is clicked.

        Parameters
        ----------
        callback : callable
            Zero-argument function.
        """
        self._reset_view_callbacks.append(callback)

    def on_request_snapshot(self, callback: Callable[[], None]):
        """Register a callback invoked when the Snapshot button is clicked.

        Parameters
        ----------
        callback : callable
            Zero-argument function.
        """
        self._snapshot_callbacks.append(callback)

    def _on_change(self):
        if self.save_config_on_exit:
            self.save_settings()
        for cb in self._display_options_callbacks:
            cb()

    def _on_reset_view_clicked(self, _btn):
        for cb in self._reset_view_callbacks:
            cb()

    def _on_snapshot_clicked(self, _btn):
        for cb in self._snapshot_callbacks:
            cb()

    def _on_settings_clicked(self, _btn):
        if self.settings_panel.layout.display == "none":
            self.settings_panel.sync_from_options(self.display_options())
            self.settings_panel.layout.display = "flex"
        else:
            self.settings_panel.layout.display = "none"

    def _on_settings_applied(self, updated: DisplayOptions):
        self._background_color = updated.background_color
        self._sites_color = updated.default_site_color
        self._exchanges_color = updated.default_exchange_color
        self._selected_color = updated.selected_color
        self._hover_color = updated.hover_color
        self._selected_hover_color = updated.selected_hover_color
        self._anisotropy_color = updated.anisotropy_color
        self._on_change()

    def display_options(self) -> DisplayOptions:
        """Build a ``DisplayOptions`` instance from the current state of the toolbar controls.

        Returns
        -------
        DisplayOptions
            Current display options.
        """
        return DisplayOptions(
            show_sites=bool(self.show_sites.value),
            show_exchanges=bool(self.show_exchanges.value),
            show_anisotropies=bool(self.show_anisotropies.value),
            show_unit_cell=bool(self.show_cell.value),
            show_supercell=bool(self.show_supercell.value),
            show_spins_not_atoms=bool(self.atom_or_spins.value),
            show_nonmagnetic_atoms=bool(self.show_nonmagnetic.value),
            use_atomic_radii=bool(self.scale_atoms.value),
            atom_spin_scaling=self.spin_scale_slider.value(),
            exchange_scaling=self.exchange_scale_slider.value(),
            anisotropy_scaling=self.anisotropy_scale_slider.value(),
            show_cartesian_axes=bool(self.show_cartesian_axes.value),
            show_lattice_axes=bool(self.show_lattice_axes.value),
            orthogonal_lattice_axes=bool(self.orthogonal_lattice_axes.value),
            perspective=bool(self.orthographic.value),
            prettify=bool(self.prettify.value),
            background_color=self._background_color,
            default_site_color=self._sites_color,
            default_exchange_color=self._exchanges_color,
            selected_color=self._selected_color,
            hover_color=self._hover_color,
            selected_hover_color=self._selected_hover_color,
            anisotropy_color=self._anisotropy_color,
        )

    def save_settings(self):
        """Serialise and write current display options to ``self.settings_filename``."""
        with open(self.settings_filename, "w", encoding="utf-8") as file:
            file.write(self.display_options().serialise())
