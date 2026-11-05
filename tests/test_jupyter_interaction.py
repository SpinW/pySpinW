"""Unit tests for Jupyter canvas interaction handler using real CrystalViewerWidget."""

from unittest.mock import MagicMock
import numpy as np
import pytest

from pyspinw.gui.jupyter.crystalview_jupyter import CrystalViewerWidget
from pyspinw.gui.jupyter.interaction import CanvasInteractionHandler
from pyspinw.gui.rendermodel import RenderModel
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.interface import generate_structure
from pyspinw.symmetry.unitcell import UnitCell


@pytest.fixture
def viewer_widget():
    """Create a real CrystalViewerWidget instance with a minimal RenderModel."""
    unit_cell = UnitCell(3.0, 3.0, 3.0)
    sites = generate_structure(
        unit_cell,
        positions=[[0.0, 0.0, 0.0]],
        spins=[[0.0, 0.0, 1.0]],
        names=["Fe"],
    )
    hamiltonian = Hamiltonian(sites, [])
    render_model = RenderModel(hamiltonian)
    return CrystalViewerWidget(render_model)


def test_orbit_rotation(viewer_widget):
    """Verify left-drag updates mouse_rotation and commits on mouseup."""
    widget = viewer_widget
    handler = widget.interaction

    render_calls = 0
    orig_render = widget.render

    def mock_render():
        nonlocal render_calls
        render_calls += 1
        orig_render()

    widget.render = mock_render

    # Mouse down
    handler.handle_event({"type": "mousedown", "button": 0, "clientX": 100, "clientY": 100})
    assert handler.is_dragging
    assert handler.button == 0

    # Mouse move (dx=20, dy=10)
    handler.handle_event({"type": "mousemove", "clientX": 120, "clientY": 110})
    assert render_calls == 1
    # Check that mouse_rotation changed from identity
    assert not np.allclose(widget.mouse_rotation, np.eye(3))

    # Mouse up
    handler.handle_event({"type": "mouseup", "shiftKey": False})
    assert not handler.is_dragging
    assert render_calls == 2
    # Committed into view_rotation
    assert not np.allclose(widget.view_rotation, np.eye(3))
    # mouse_rotation reset to identity
    assert np.allclose(widget.mouse_rotation, np.eye(3))


def test_pan_movement(viewer_widget):
    """Verify right-drag updates mouse_origin and commits on mouseup."""
    widget = viewer_widget
    handler = widget.interaction
    initial_origin = widget.view_origin.copy()

    render_calls = 0
    orig_render = widget.render

    def mock_render():
        nonlocal render_calls
        render_calls += 1
        orig_render()

    widget.render = mock_render

    # Right mouse down (button 2)
    handler.handle_event({"type": "mousedown", "button": 2, "clientX": 50, "clientY": 50})
    assert handler.button == 2

    # Mouse move (dx=30, dy=40)
    handler.handle_event({"type": "mousemove", "clientX": 80, "clientY": 90})
    assert render_calls == 1
    assert not np.allclose(widget.mouse_origin, np.zeros(3))

    # Mouse up
    handler.handle_event({"type": "mouseup"})
    assert render_calls == 2
    assert not np.allclose(widget.view_origin, initial_origin)
    assert np.allclose(widget.mouse_origin, np.zeros(3))


def test_click_detection(viewer_widget):
    """Verify that a click without significant drag triggers on_canvas_click."""
    widget = viewer_widget
    handler = widget.interaction
    widget.on_canvas_click = MagicMock()

    handler.handle_event({"type": "mousedown", "button": 0, "clientX": 100, "clientY": 100})
    # Tiny move (within dx^2 + dy^2 < 100)
    handler.handle_event({"type": "mousemove", "clientX": 102, "clientY": 101})
    handler.handle_event({"type": "mouseup", "shiftKey": True})

    widget.on_canvas_click.assert_called_once_with(shift=True)


def test_wheel_zoom(viewer_widget):
    """Verify wheel zoom changes view_radius and respects min_view_radius."""
    widget = viewer_widget
    handler = CanvasInteractionHandler(widget, min_view_radius=0.5)

    initial_radius = widget.view_radius
    # Scroll down / zoom out (positive deltaY)
    handler.handle_event({"type": "wheel", "deltaY": 100.0})
    assert widget.view_radius > initial_radius

    # Scroll up / zoom in (negative deltaY)
    handler.handle_event({"type": "wheel", "deltaY": -200.0})
    assert widget.view_radius < initial_radius

    # Zoom in excessively: should clamp to min_view_radius
    handler.handle_event({"type": "wheel", "deltaY": -10000.0})
    assert np.isclose(widget.view_radius, 0.5)
