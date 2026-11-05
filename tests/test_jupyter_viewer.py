"""Unit tests for Jupyter notebook composite viewer, toolbar controls, and view_jupyter."""

from unittest.mock import MagicMock
import numpy as np
import pytest

from pyspinw.gui.displayoptions import DisplayOptions
from pyspinw.gui.jupyter import (
    CrystalViewerWidget,
    DisplayControlsWidget,
    ViewerJupyter,
    view_jupyter,
)
from pyspinw.gui.rendermodel import RenderModel
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.interface import generate_structure
from pyspinw.structure import Structure
from pyspinw.symmetry.unitcell import UnitCell


@pytest.fixture
def sample_hamiltonian():
    """Create a minimal Hamiltonian with sites and unit cell for testing."""
    unit_cell = UnitCell(3.0, 4.0, 5.0, 90.0, 90.0, 90.0)
    sites = generate_structure(
        unit_cell,
        positions=[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
        spins=[[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]],
        names=["Fe1", "Fe2"],
    )
    return Hamiltonian(sites, [])


@pytest.fixture
def viewer_widget(sample_hamiltonian):
    """Create a CrystalViewerWidget instance."""
    render_model = RenderModel(sample_hamiltonian)
    return CrystalViewerWidget(render_model)


def test_align_view_axes(viewer_widget):
    """Test align_view with various crystallographic and Cartesian axes."""
    # Test Cartesian alignments
    viewer_widget.align_view("x")
    assert not np.allclose(viewer_widget.view_rotation, np.eye(3))

    viewer_widget.align_view("z")
    assert np.allclose(viewer_widget.view_rotation, np.eye(3))

    # Test lattice alignments
    viewer_widget.align_view("a")
    assert viewer_widget.view_rotation.shape == (3, 3)

    viewer_widget.align_view("b")
    assert viewer_widget.view_rotation.shape == (3, 3)

    viewer_widget.align_view("c")
    assert viewer_widget.view_rotation.shape == (3, 3)

    viewer_widget.align_view("a*")
    assert viewer_widget.view_rotation.shape == (3, 3)

    # Test invalid axis
    with pytest.raises(ValueError, match="Unknown axis"):
        viewer_widget.align_view("invalid_axis")


def test_display_controls_widget(viewer_widget):
    """Test toolbar toggles and sliders sync with viewer.display_options."""
    controls = DisplayControlsWidget(viewer=viewer_widget)

    render_calls = 0
    orig_render = viewer_widget.render

    def count_render():
        nonlocal render_calls
        render_calls += 1
        orig_render()

    viewer_widget.render = count_render

    # Test toggles
    controls.btn_sites.value = False
    assert viewer_widget.display_options.show_sites is False
    assert render_calls == 1

    controls.btn_spins.value = False
    assert viewer_widget.display_options.show_spins_not_atoms is False
    assert render_calls == 2

    controls.btn_exchanges.value = False
    assert viewer_widget.display_options.show_exchanges is False
    assert render_calls == 3

    # Test sliders
    controls.slider_spin_scale.value = 2.5
    assert np.isclose(viewer_widget.display_options.atom_spin_scaling, 2.5)
    assert render_calls == 4

    controls.slider_exchange_scale.value = 1.8
    assert np.isclose(viewer_widget.display_options.exchange_scaling, 1.8)
    assert render_calls == 5

    # Test reset button
    viewer_widget.view_radius = 50.0
    controls.btn_reset.click()
    assert np.isclose(viewer_widget.view_radius, 10.0)

    # Test align dropdown
    controls.align_dropdown.value = "c"
    assert controls.align_dropdown.value == "Align view..."


def test_viewer_jupyter_composite(sample_hamiltonian):
    """Test composite ViewerJupyter initialization and structure info banner."""
    viewer_jupyter = ViewerJupyter(sample_hamiltonian, width=600, height=400)

    assert isinstance(viewer_jupyter.viewer, CrystalViewerWidget)
    assert isinstance(viewer_jupyter.controls, DisplayControlsWidget)
    assert len(viewer_jupyter.children) == 3

    # Info label contains lattice and site counts
    info_html = viewer_jupyter.info_label.value
    assert "3.000" in info_html
    assert "4.000" in info_html
    assert "5.000" in info_html
    assert "Sites:" in info_html


def test_view_jupyter_convenience(sample_hamiltonian):
    """Test view_jupyter helper with Hamiltonian, Structure, and invalid input."""
    # 1. With Hamiltonian
    v1 = view_jupyter(sample_hamiltonian)
    assert isinstance(v1, ViewerJupyter)

    # 2. With Structure
    structure = sample_hamiltonian.structure
    assert isinstance(structure, Structure)
    v2 = view_jupyter(structure)
    assert isinstance(v2, ViewerJupyter)

    # 3. With invalid input
    with pytest.raises(TypeError, match="view_jupyter expected Hamiltonian or Structure"):
        view_jupyter("invalid_object")  # type: ignore

