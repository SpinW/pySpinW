"""Interactive controls and toolbar for Jupyter Notebook crystal visualization."""

from __future__ import annotations

from typing import TYPE_CHECKING, Callable
import ipywidgets as widgets

from pyspinw.gui.displayoptions import DisplayOptions

if TYPE_CHECKING:
    from pyspinw.gui.jupyter.crystalview_jupyter import CrystalViewerWidget


class DisplayControlsWidget(widgets.VBox):
    """Jupyter toolbar widget for interacting with DisplayOptions and camera alignment.

    Parameters
    ----------
    viewer : CrystalViewerWidget or None, optional
        The crystal viewer widget to control.
    on_change : callable or None, optional
        Optional callback called whenever display settings change.
    """

    def __init__(
        self,
        viewer: CrystalViewerWidget | None = None,
        on_change: Callable[[], None] | None = None,
    ):
        super().__init__()
        self.viewer = viewer
        self.on_change = on_change

        opts = viewer.display_options if viewer is not None else DisplayOptions()

        # 1. Visibility Toggles
        button_layout = widgets.Layout(width="auto", height="32px", margin="2px 4px 2px 0px")

        self.btn_sites = widgets.ToggleButton(
            value=opts.show_sites,
            description="Sites",
            tooltip="Show atomic/spin sites",
            layout=button_layout,
        )
        self.btn_spins = widgets.ToggleButton(
            value=opts.show_spins_not_atoms,
            description="Spins",
            tooltip="Toggle between showing spin vectors and atomic spheres",
            layout=button_layout,
        )
        self.btn_exchanges = widgets.ToggleButton(
            value=opts.show_exchanges,
            description="Bonds",
            tooltip="Show magnetic exchange bonds",
            layout=button_layout,
        )
        self.btn_anisotropies = widgets.ToggleButton(
            value=opts.show_anisotropies,
            description="Aniso",
            tooltip="Show single-ion anisotropy ellipsoids",
            layout=button_layout,
        )
        self.btn_cell = widgets.ToggleButton(
            value=opts.show_unit_cell,
            description="Unit Cell",
            tooltip="Show unit cell wireframe",
            layout=button_layout,
        )
        self.btn_supercell = widgets.ToggleButton(
            value=opts.show_supercell,
            description="Supercell",
            tooltip="Show supercell wireframe",
            layout=button_layout,
        )
        self.btn_axes = widgets.ToggleButton(
            value=opts.show_cartesian_axes or opts.show_lattice_axes,
            description="Axes",
            tooltip="Show corner coordinate axes",
            layout=button_layout,
        )

        toggles_row = widgets.HBox(
            [
                self.btn_sites,
                self.btn_spins,
                self.btn_exchanges,
                self.btn_anisotropies,
                self.btn_cell,
                self.btn_supercell,
                self.btn_axes,
            ],
            layout=widgets.Layout(flex_flow="row wrap", margin="2px 0px"),
        )

        # 2. Camera Controls & Alignment
        self.btn_reset = widgets.Button(
            description="Reset View",
            icon="undo",
            tooltip="Reset camera view to default orientation",
            layout=widgets.Layout(width="110px", height="32px", margin="2px 4px 2px 0px"),
        )
        self.btn_perspective = widgets.ToggleButton(
            value=opts.perspective,
            description="Perspective",
            tooltip="Toggle perspective vs orthographic projection",
            layout=button_layout,
        )
        self.align_dropdown = widgets.Dropdown(
            options=["Align view...", "a", "b", "c", "a*", "b*", "c*", "x", "y", "z"],
            value="Align view...",
            description="",
            layout=widgets.Layout(width="130px", height="32px", margin="2px 4px 2px 0px"),
        )

        camera_row = widgets.HBox(
            [self.btn_reset, self.btn_perspective, self.align_dropdown],
            layout=widgets.Layout(flex_flow="row wrap", margin="2px 0px"),
        )

        # 3. Scaling Sliders
        slider_layout = widgets.Layout(width="220px", margin="2px 10px 2px 0px")
        self.slider_spin_scale = widgets.FloatSlider(
            value=opts.atom_spin_scaling,
            min=0.1,
            max=5.0,
            step=0.1,
            description="Spin:",
            continuous_update=False,
            layout=slider_layout,
        )
        self.slider_exchange_scale = widgets.FloatSlider(
            value=opts.exchange_scaling,
            min=0.1,
            max=5.0,
            step=0.1,
            description="Bond:",
            continuous_update=False,
            layout=slider_layout,
        )
        self.slider_aniso_scale = widgets.FloatSlider(
            value=opts.anisotropy_scaling,
            min=0.1,
            max=5.0,
            step=0.1,
            description="Aniso:",
            continuous_update=False,
            layout=slider_layout,
        )

        sliders_row = widgets.HBox(
            [self.slider_spin_scale, self.slider_exchange_scale, self.slider_aniso_scale],
            layout=widgets.Layout(flex_flow="row wrap", margin="2px 0px"),
        )

        # Wire up event observers
        self.btn_sites.observe(self._on_toggle_change, names="value")
        self.btn_spins.observe(self._on_toggle_change, names="value")
        self.btn_exchanges.observe(self._on_toggle_change, names="value")
        self.btn_anisotropies.observe(self._on_toggle_change, names="value")
        self.btn_cell.observe(self._on_toggle_change, names="value")
        self.btn_supercell.observe(self._on_toggle_change, names="value")
        self.btn_axes.observe(self._on_toggle_change, names="value")
        self.btn_perspective.observe(self._on_toggle_change, names="value")

        self.slider_spin_scale.observe(self._on_slider_change, names="value")
        self.slider_exchange_scale.observe(self._on_slider_change, names="value")
        self.slider_aniso_scale.observe(self._on_slider_change, names="value")

        self.btn_reset.on_click(self._on_reset_clicked)
        self.align_dropdown.observe(self._on_align_selected, names="value")

        self.children = [toggles_row, camera_row, sliders_row]

    def _sync_to_viewer(self) -> None:
        """Apply current control values to viewer display_options and request redraw."""
        if self.viewer is None:
            return

        opts = self.viewer.display_options
        opts.show_sites = self.btn_sites.value
        opts.show_spins_not_atoms = self.btn_spins.value
        opts.show_exchanges = self.btn_exchanges.value
        opts.show_anisotropies = self.btn_anisotropies.value
        opts.show_unit_cell = self.btn_cell.value
        opts.show_supercell = self.btn_supercell.value
        opts.show_cartesian_axes = self.btn_axes.value
        opts.show_lattice_axes = self.btn_axes.value
        opts.perspective = self.btn_perspective.value

        opts.atom_spin_scaling = float(self.slider_spin_scale.value)
        opts.exchange_scaling = float(self.slider_exchange_scale.value)
        opts.anisotropy_scaling = float(self.slider_aniso_scale.value)

        self.viewer.request_render()

        if self.on_change is not None:
            self.on_change()

    def _on_toggle_change(self, change) -> None:
        self._sync_to_viewer()

    def _on_slider_change(self, change) -> None:
        self._sync_to_viewer()

    def _on_reset_clicked(self, btn) -> None:
        if self.viewer is not None:
            self.viewer.reset_view()

    def _on_align_selected(self, change) -> None:
        val = change["new"]
        if val and val != "Align view..." and self.viewer is not None:
            self.viewer.align_view(val)
            self.align_dropdown.value = "Align view..."

