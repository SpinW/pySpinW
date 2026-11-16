"""Tests for the Jupyter ipywidgets / ipywebgl viewer components."""

import os
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from pyspinw import (
    AxisMagnitudeAnisotropy,
    DisplayOptions,
    Hamiltonian,
    HeisenbergExchange,
    Site,
    Structure,
    TiledSupercell,
    UnitCell,
    generate_exchanges,
    view,
)
from pyspinw.gui.jupyter import (
    JupyterCrystalView,
    JupyterSettingsPanel,
    JupyterTextDisplay,
    JupyterToolbar,
    JupyterViewer,
)
from pyspinw.gui.jupyter.toolbar import hex_to_rgb, rgb_to_hex
from pyspinw.gui.rendering.models.arrow import Arrow
from pyspinw.gui.rendering.models.sphere import Sphere
from pyspinw.gui.rendering.models.tube import Tube
from pyspinw.gui.rendering.models.wrireframe_cube import WireframeCube


@pytest.fixture
def sample_hamiltonian() -> Hamiltonian:
    """Create a sample Hamiltonian with magnetic/non-magnetic sites, exchanges, and anisotropy."""
    unit_cell = UnitCell(1, 1, 1, gamma=60)
    x = Site(0, 0, 0, 0, 0, 1, name="X")
    y = Site(0.5, 0, 0, 0, 0, 1, name="Y")
    z = Site(0, 0.5, 0, 0, 0, 1, name="Z")
    w = Site(0.5, 0.5, 0, 0, 0, 0, name="W")

    structure = Structure(
        [x, y, z, w],
        unit_cell=unit_cell,
        supercell=TiledSupercell(scaling=(2, 2, 1)),
    )
    exchanges = generate_exchanges(
        sites=[x, y, z],
        unit_cell=unit_cell,
        max_distance=0.6,
        exchange_type=HeisenbergExchange,
        j=-1.0,
    )
    anisotropies = [AxisMagnitudeAnisotropy(x, -0.5, [0, 0, 1])]
    return Hamiltonian(structure, exchanges, anisotropies)


def test_headless_geometry_builders():
    """Verify that mesh geometry builders run purely in NumPy without an OpenGL context."""
    arrow_sections = Arrow.build_geometry()
    assert len(arrow_sections) == 4
    for section in arrow_sections:
        assert section.ndim == 2
        assert section.shape[1] == 6
        assert section.dtype == np.float32

    sphere_data = Sphere.build_geometry(3, radius=0.5, padding=0.05)
    assert sphere_data.ndim == 2
    assert sphere_data.shape[1] == 3
    assert sphere_data.dtype == np.float32

    tube_data = Tube.build_geometry(n_points=21, padding=0.1)
    assert tube_data.shape == (20 * 6, 6)
    assert tube_data.dtype == np.float32

    cube_data = WireframeCube.build_geometry()
    assert cube_data.shape == (24, 3)
    assert cube_data.dtype == np.float32


def test_color_conversion_helpers():
    """Verify conversion between RGB float triples and hex color strings."""
    hex_str = rgb_to_hex((1.0, 0.5, 0.0))
    assert hex_str == "#ff8000"
    rgb = hex_to_rgb(hex_str)
    assert np.allclose(rgb, (1.0, 128 / 255.0, 0.0), atol=1e-3)


def test_jupyter_viewer_parameter_validation(sample_hamiltonian):
    """Verify constructor validation on JupyterViewer matches desktop Viewer."""
    with pytest.raises(ValueError, match="Expected `copies` to be a triple of ints"):
        JupyterViewer(sample_hamiltonian, copies=(1, 2), save_config_on_exit=False)

    with pytest.raises(ValueError, match="Expected `rotation_supercell_expansion_max` to be a triple of ints"):
        JupyterViewer(sample_hamiltonian, rotation_supercell_expansion_max=(1, 2), save_config_on_exit=False)

    with pytest.raises(ValueError, match="Expected initial_rotation to be orthogonal"):
        JupyterViewer(sample_hamiltonian, initial_rotation=2.0 * np.eye(3), save_config_on_exit=False)

    with pytest.raises(ValueError, match="Expected a rotation, not a rotoreflection"):
        JupyterViewer(sample_hamiltonian, initial_rotation=-np.eye(3), save_config_on_exit=False)

    with pytest.raises(ValueError, match="Expected initial distance to be strictly positive"):
        JupyterViewer(sample_hamiltonian, initial_distance=0.0, save_config_on_exit=False)


def test_jupyter_viewer_toolbar_and_selection_sync(sample_hamiltonian, tmp_path):
    """Verify toolbar updates, settings panel, text selection, hover, and snapshot saving."""
    viewer = JupyterViewer(sample_hamiltonian, save_config_on_exit=False)

    assert isinstance(viewer.toolbar, JupyterToolbar)
    assert isinstance(viewer.viewer, JupyterCrystalView)
    assert isinstance(viewer.text_display, JupyterTextDisplay)
    assert isinstance(viewer.toolbar.settings_panel, JupyterSettingsPanel)

    # Toggle display options on the toolbar and verify synchronization with JupyterCrystalView
    viewer.toolbar.show_sites.value = False
    assert viewer.viewer.display_options.show_sites is False

    viewer.toolbar.orthographic.value = False
    assert viewer.viewer.display_options.perspective is False

    viewer.toolbar.show_lattice_axes.value = True
    viewer.toolbar.orthogonal_lattice_axes.value = True
    assert viewer.viewer.display_options.show_lattice_axes is True
    assert viewer.viewer.display_options.orthogonal_lattice_axes is True

    viewer.toolbar.spin_scale_slider.slider.value = 0.75
    assert np.isclose(viewer.viewer.display_options.atom_spin_scaling, 0.75)

    # Settings panel color update
    viewer.toolbar._on_settings_clicked(None)
    assert viewer.toolbar.settings_panel.layout.display == "flex"
    viewer.toolbar.settings_panel.background_picker.value = "#112233"
    viewer.toolbar.settings_panel._on_ok_clicked(None)
    assert viewer.toolbar.settings_panel.layout.display == "none"
    assert np.allclose(viewer.viewer.display_options.background_color, hex_to_rgb("#112233"))

    # Select sites via text_display site_tree and verify 3D viewer selection updates
    first_site_node = viewer.text_display.site_tree.nodes[0]
    first_site_node.selected = True
    assert len(viewer.viewer.current_selection) > 0
    assert viewer.viewer.current_selection == first_site_node.ids
    assert viewer.viewer.current_selection == viewer.text_display.current_selection

    # Deselect site node
    first_site_node.selected = False
    assert viewer.viewer.current_selection == []

    # Click on element from 3D viewer side and verify text_display updates and expands parent
    viewer.viewer.click_on_element(1, shift=False)
    assert viewer.text_display.current_selection == [1]
    assert viewer.text_display.site_render_id_to_row_item[1][0].opened is True
    assert viewer.text_display.site_render_id_to_row_item[1][1].selected is True
    viewer.viewer.click_on_element(2, shift=True)
    assert viewer.text_display.current_selection == [1, 2]
    viewer.viewer.click_on_element(0, shift=False)
    assert viewer.text_display.current_selection == []
    assert viewer.text_display.site_tree.selected_nodes == ()

    # Hover synchronization and row bolding
    viewer.text_display.set_hover([1])
    assert "font-weight:bold" in viewer.text_display.site_render_id_to_row_item[1][1].name
    viewer.text_display.site_tree.hover_ids = [1]
    viewer.text_display.on_hover_changed()
    assert viewer.viewer.hover_ids == [1]

    # Reset view
    viewer.toolbar._on_reset_view_clicked(None)
    assert np.isclose(viewer.viewer.view_radius, 10.0)

    # Snapshot saving
    snap_path = tmp_path / "test_snap.png"
    viewer.on_snapshot_requested()
    assert viewer.snapshot_bar.layout.display == "flex"
    viewer.snapshot_filename_input.value = str(snap_path)
    viewer._on_snapshot_confirm(None)
    assert viewer.snapshot_bar.layout.display == "none"
    assert os.path.exists(snap_path)


def test_view_dispatch_jupyter_mode(sample_hamiltonian):
    """Verify pyspinw.view (show_object) returns a JupyterViewer when jupyter=True."""
    with patch("IPython.display.display", MagicMock()) as mock_display:
        widget = view(sample_hamiltonian.structure, jupyter=True)
        assert isinstance(widget, JupyterViewer)
        mock_display.assert_called_once_with(widget)
