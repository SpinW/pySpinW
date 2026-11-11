"""Main Jupyter notebook viewer widget combining toolbar, WebGL 3D canvas, and text tables."""

from imageio import imwrite
import ipywidgets as widgets
from IPython.display import display
import numpy as np

from pyspinw.gui.displayoptions import DisplayOptions
from pyspinw.gui.jupyter.crystalview import JupyterCrystalView
from pyspinw.gui.jupyter.textdisplay import JupyterTextDisplay
from pyspinw.gui.jupyter.toolbar import JupyterToolbar
from pyspinw.gui.rendermodel import RenderModel
from pyspinw.hamiltonian import Hamiltonian

_unique_id_counter = -1


def _generate_unique_id() -> int:
    """Generate a unique integer ID for each Jupyter viewer instance.

    Returns
    -------
    int
        Unique integer identifier.
    """
    global _unique_id_counter  # noqa: PLW0603
    _unique_id_counter += 1
    return _unique_id_counter


class JupyterViewer(widgets.VBox):
    """Main Jupyter notebook viewer widget using ``ipywidgets`` and ``ipywebgl``.

    Parameters
    ----------
    hamiltonian : Hamiltonian
        Hamiltonian to display in the viewer.
    copies : tuple of (int, int, int), optional
        Supercell scaling along each lattice vector. Default is ``(1, 1, 1)``.
    rotation_supercell_expansion_max : tuple of (int, int, int), optional
        Maximum denominator along each axis when approximating a
        ``RotationSupercell`` propagation vector as a commensurate cell.
        Default is ``(10, 10, 10)``.
    propagation_vector_commensurability_tolerance : float, optional
        Tolerance for matching a propagation vector component to a rational
        fraction. Default is ``1e-4``.
    parent : object or None, optional
        Unused parent reference kept for API compatibility with ``Viewer``.
        Default is ``None``.
    initial_rotation : (3, 3) numpy.ndarray or None, optional
        Initial 3x3 view rotation matrix. Default is ``None``.
    initial_distance : float or None, optional
        Initial camera distance from the view origin. Default is ``None``.
    display_options : DisplayOptions or None, optional
        Initial display options. Default is ``None``.
    save_config_on_exit : bool, optional
        Whether to save display configuration when display options change.
        Default is ``True``.

    Raises
    ------
    ValueError
        If ``copies`` or ``rotation_supercell_expansion_max`` is not a triple
        of ints, if ``initial_rotation`` is not a valid rotation matrix, or if
        ``initial_distance`` is not strictly positive.
    """

    def __init__(
        self,
        hamiltonian: Hamiltonian,
        copies: tuple[int, int, int] = (1, 1, 1),
        rotation_supercell_expansion_max: tuple[int, int, int] = (10, 10, 10),
        propagation_vector_commensurability_tolerance: float = 1e-4,
        parent=None,
        initial_rotation: np.ndarray | None = None,
        initial_distance: float | None = None,
        display_options: DisplayOptions | None = None,
        save_config_on_exit: bool = True,
    ):
        if len(copies) != 3 or not all(isinstance(x, int) for x in copies):
            raise ValueError("Expected `copies` to be a triple of ints")

        if len(rotation_supercell_expansion_max) != 3 or not all(
            isinstance(x, int) for x in rotation_supercell_expansion_max
        ):
            raise ValueError("Expected `rotation_supercell_expansion_max` to be a triple of ints")

        if initial_rotation is not None:
            if not np.allclose(initial_rotation @ initial_rotation.T, np.eye(3)):
                raise ValueError("Expected initial_rotation to be orthogonal")

            if not np.isclose(np.linalg.det(initial_rotation), 1):
                raise ValueError("Expected a rotation, not a rotoreflection")

        if initial_distance is not None and initial_distance <= 0:
            raise ValueError("Expected initial distance to be strictly positive")

        self.save_config_on_exit = save_config_on_exit
        self._unique_id = _generate_unique_id()

        self.render_model = RenderModel(
            hamiltonian,
            copies,
            rotation_supercell_expansion_max,
            propagation_vector_commensurability_tolerance,
        )

        self.toolbar = JupyterToolbar(
            initial_display_options=display_options,
            save_config_on_exit=self.save_config_on_exit,
            max_size_object_scaling=float(np.min(hamiltonian.structure.unit_cell.abc)),
        )

        self.viewer = JupyterCrystalView(self.render_model, initial_rotation, initial_distance)
        self.text_display = JupyterTextDisplay(self.render_model)

        # Inline snapshot save controls for notebook environments
        self.snapshot_filename_input = widgets.Text(
            value="snapshot.png",
            description="Save to:",
            style={"description_width": "65px"},
            layout=widgets.Layout(width="260px"),
        )
        self.snapshot_save_button = widgets.Button(
            description="Save",
            button_style="primary",
            layout=widgets.Layout(width="70px"),
        )
        self.snapshot_cancel_button = widgets.Button(
            description="Cancel",
            layout=widgets.Layout(width="70px"),
        )
        self.snapshot_bar = widgets.HBox(
            [self.snapshot_filename_input, self.snapshot_save_button, self.snapshot_cancel_button],
            layout=widgets.Layout(display="none", align_items="center", margin="2px 0px"),
        )
        self.snapshot_save_button.on_click(self._on_snapshot_confirm)
        self.snapshot_cancel_button.on_click(self._on_snapshot_cancel)

        self.split_pane = widgets.HBox(
            [self.viewer.widget, self.text_display],
            layout=widgets.Layout(width="100%", align_items="flex-start"),
        )

        super().__init__(
            children=[self.toolbar, self.snapshot_bar, self.split_pane],
            layout=widgets.Layout(width="100%"),
        )

        # Wire up toolbar settings
        self.viewer.display_options = self.toolbar.display_options()
        self.toolbar.on_display_options_changed(self.on_display_options_changed)
        self.toolbar.on_request_view_reset(self.on_reset_view_requested)
        self.toolbar.on_request_snapshot(self.on_snapshot_requested)

        # Wire up text and graphics selection / hover
        self.text_display.on_hover_changed(self.on_text_hover_changed)
        self.text_display.on_selection_changed(self.on_text_selection_changed)
        self.viewer.on_hover_changed(self.on_render_hover_changed)
        self.viewer.on_selection_changed(self.on_render_selection_changed)

    def on_snapshot_requested(self):
        """Toggle the inline snapshot filename bar when the snapshot button is clicked."""
        if self.snapshot_bar.layout.display == "none":
            self.snapshot_bar.layout.display = "flex"
        else:
            self.snapshot_bar.layout.display = "none"

    def _on_snapshot_confirm(self, _btn):
        filename = self.snapshot_filename_input.value.strip()
        if filename:
            self.save_snapshot(filename)
        self.snapshot_bar.layout.display = "none"

    def _on_snapshot_cancel(self, _btn):
        self.snapshot_bar.layout.display = "none"

    def save_snapshot(self, filename: str):
        """Write the current viewer image data to a file.

        Parameters
        ----------
        filename : str
            Output image file path.
        """
        data = self.viewer.snapshot()
        imwrite(filename, data)

    def on_display_options_changed(self):
        """Update the 3D viewer when toolbar display options change."""
        self.viewer.display_options = self.toolbar.display_options()

    def on_text_hover_changed(self):
        """Update the 3D viewer when hovered text items change."""
        self.viewer.hover_ids = self.text_display.hover_ids

    def on_text_selection_changed(self):
        """Update the 3D viewer when selected text items change."""
        self.viewer.current_selection = self.text_display.current_selection

    def on_render_hover_changed(self):
        """Update the text display when hovered 3D items change."""
        self.text_display.set_hover(self.viewer.hover_ids)

    def on_render_selection_changed(self):
        """Update the text display when selected 3D items change."""
        self.text_display.set_selection(self.viewer.current_selection)

    def on_reset_view_requested(self):
        """Reset the 3D camera view when requested from the toolbar."""
        self.viewer.reset_view()
