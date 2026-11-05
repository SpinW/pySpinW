"""Unit tests for Jupyter canvas interaction handler."""

import numpy as np
import pytest

from pyspinw.gui.jupyter.interaction import CanvasInteractionHandler


class DummyViewerWidget:
    """Mock widget holding view parameters for interaction testing."""

    def __init__(self):
        self.view_origin = np.array([0.0, 0.0, 0.0], dtype=float)
        self.view_rotation = np.eye(3, dtype=float)
        self.view_radius = 10.0
        self.mouse_rotation = np.eye(3, dtype=float)
        self.mouse_origin = np.zeros(3, dtype=float)

        self.render_count = 0
        self.clicked = False
        self.clicked_shift = False

    def request_render(self):
        self.render_count += 1

    def on_canvas_click(self, shift=False):
        self.clicked = True
        self.clicked_shift = shift


def test_orbit_rotation():
    """Verify left-drag updates mouse_rotation and commits on mouseup."""
    widget = DummyViewerWidget()
    handler = CanvasInteractionHandler(widget)

    # Mouse down
    handler.handle_event({"type": "mousedown", "button": 0, "clientX": 100, "clientY": 100})
    assert handler.is_dragging
    assert handler.button == 0

    # Mouse move (dx=20, dy=10)
    handler.handle_event({"type": "mousemove", "clientX": 120, "clientY": 110})
    assert widget.render_count == 1
    # Check that mouse_rotation changed from identity
    assert not np.allclose(widget.mouse_rotation, np.eye(3))

    # Mouse up
    handler.handle_event({"type": "mouseup", "shiftKey": False})
    assert not handler.is_dragging
    assert widget.render_count == 2
    # Committed into view_rotation
    assert not np.allclose(widget.view_rotation, np.eye(3))
    # mouse_rotation reset to identity
    assert np.allclose(widget.mouse_rotation, np.eye(3))


def test_pan_movement():
    """Verify right-drag updates mouse_origin and commits on mouseup."""
    widget = DummyViewerWidget()
    handler = CanvasInteractionHandler(widget)

    # Right mouse down (button 2)
    handler.handle_event({"type": "mousedown", "button": 2, "clientX": 50, "clientY": 50})
    assert handler.button == 2

    # Mouse move (dx=30, dy=40)
    handler.handle_event({"type": "mousemove", "clientX": 80, "clientY": 90})
    assert widget.render_count == 1
    assert not np.allclose(widget.mouse_origin, np.zeros(3))

    # Mouse up
    handler.handle_event({"type": "mouseup"})
    assert widget.render_count == 2
    assert not np.allclose(widget.view_origin, np.zeros(3))
    assert np.allclose(widget.mouse_origin, np.zeros(3))


def test_click_detection():
    """Verify that a click without significant drag triggers on_canvas_click."""
    widget = DummyViewerWidget()
    handler = CanvasInteractionHandler(widget)

    handler.handle_event({"type": "mousedown", "button": 0, "clientX": 100, "clientY": 100})
    # Tiny move (within dx^2 + dy^2 < 100)
    handler.handle_event({"type": "mousemove", "clientX": 102, "clientY": 101})
    handler.handle_event({"type": "mouseup", "shiftKey": True})

    assert widget.clicked
    assert widget.clicked_shift


def test_wheel_zoom():
    """Verify wheel zoom changes view_radius and respects min_view_radius."""
    widget = DummyViewerWidget()
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

