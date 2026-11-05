"""Canvas mouse and interaction handler for Jupyter notebooks using ipyevents."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable
import numpy as np

from pyspinw.util import rotation_matrix

if TYPE_CHECKING:
    from pyspinw.gui.jupyter.crystalview_jupyter import CrystalViewerWidget


class CanvasInteractionHandler:
    """Handles mouse orbit rotation, panning, and scroll zooming via ipyevents.

    Parameters
    ----------
    viewer_widget : CrystalViewerWidget
        The parent CrystalViewerWidget holding view parameters and rendering methods.
    mouse_angle_sensitivity : float, optional
        Sensitivity for mouse orbiting (left-click drag). Default is 0.01.
    mouse_move_sensitivity : float, optional
        Sensitivity for mouse panning (right-click drag). Default is 0.005.
    mouse_zoom_sensitivity : float, optional
        Sensitivity for wheel scroll zooming. Default is 0.001.
    min_view_radius : float, optional
        Minimum allowable camera distance. Default is 0.01.
    """

    def __init__(
        self,
        viewer_widget: Any,
        mouse_angle_sensitivity: float = 0.01,
        mouse_move_sensitivity: float = 0.005,
        mouse_zoom_sensitivity: float = 0.001,
        min_view_radius: float = 0.01,
    ):
        self.viewer_widget = viewer_widget
        self.mouse_angle_sensitivity = mouse_angle_sensitivity
        self.mouse_move_sensitivity = mouse_move_sensitivity
        self.mouse_zoom_sensitivity = mouse_zoom_sensitivity
        self.min_view_radius = min_view_radius

        self.is_dragging = False
        self.button = 0
        self.drag_start_x = 0.0
        self.drag_start_y = 0.0
        self.current_dx = 0.0
        self.current_dy = 0.0

        self.event_listener: Any = None
        self._attach_listener()

    def _attach_listener(self) -> None:
        """Attach ipyevents.Event listener to the target WebGL viewer widget."""
        try:
            from ipyevents import Event

            source = getattr(self.viewer_widget, "gl_viewer", self.viewer_widget)
            self.event_listener = Event(
                source=source,
                watched_events=["mousedown", "mousemove", "mouseup", "wheel", "contextmenu"],
                prevent_default_action=["contextmenu", "wheel"],
            )
            self.event_listener.on_dom_event(self.handle_event)
        except ImportError:
            # ipyevents not installed yet or running in headless testing
            self.event_listener = None

    def handle_event(self, event: dict[str, Any]) -> None:
        """Process a DOM event received from ipyevents.

        Parameters
        ----------
        event : dict
            Dictionary describing the DOM event (type, button, clientX, clientY, deltaY, etc.).
        """
        event_type = event.get("type", "")

        if event_type == "mousedown":
            self.on_mouse_down(
                x=float(event.get("clientX", event.get("offsetX", 0.0))),
                y=float(event.get("clientY", event.get("offsetY", 0.0))),
                button=int(event.get("button", 0)),
            )

        elif event_type == "mousemove":
            self.on_mouse_move(
                x=float(event.get("clientX", event.get("offsetX", 0.0))),
                y=float(event.get("clientY", event.get("offsetY", 0.0))),
            )

        elif event_type == "mouseup":
            shift = bool(event.get("shiftKey", False))
            self.on_mouse_up(shift=shift)

        elif event_type == "wheel":
            delta_y = float(event.get("deltaY", 0.0))
            self.on_wheel(delta_y)

        elif event_type == "contextmenu":
            # Suppress browser context menu on right-click
            pass

    def on_mouse_down(self, x: float, y: float, button: int = 0) -> None:
        """Handle mouse button down."""
        self.is_dragging = True
        self.button = button
        self.drag_start_x = x
        self.drag_start_y = y
        self.current_dx = 0.0
        self.current_dy = 0.0

        self.viewer_widget.mouse_rotation = np.eye(3, dtype=float)
        self.viewer_widget.mouse_origin = np.zeros(3, dtype=float)

    def on_mouse_move(self, x: float, y: float) -> None:
        """Handle mouse movement during drag."""
        if not self.is_dragging:
            return

        dx = x - self.drag_start_x
        dy = y - self.drag_start_y
        self.current_dx = dx
        self.current_dy = dy

        if self.button == 0:
            # Left button: Orbit rotation
            # Rotate around Y by -sensitivity*dx, and around X by -sensitivity*dy
            r_y = rotation_matrix(-self.mouse_angle_sensitivity * dx, (0.0, 1.0, 0.0))
            r_x = rotation_matrix(-self.mouse_angle_sensitivity * dy, (1.0, 0.0, 0.0))
            self.viewer_widget.mouse_rotation = r_y @ r_x
            self.viewer_widget.request_render()

        elif self.button == 2:
            # Right button: Panning view origin
            r = self.viewer_widget.view_rotation
            self.viewer_widget.mouse_origin = r @ np.array(
                [
                    -self.mouse_move_sensitivity * dx,
                    self.mouse_move_sensitivity * dy,
                    0.0,
                ],
                dtype=float,
            )
            self.viewer_widget.request_render()

    def on_mouse_up(self, shift: bool = False) -> None:
        """Handle mouse release."""
        if not self.is_dragging:
            return

        # Commit current transformations to persistent view state
        self.viewer_widget.view_rotation = (
            self.viewer_widget.view_rotation @ self.viewer_widget.mouse_rotation
        )
        self.viewer_widget.view_origin = (
            self.viewer_widget.view_origin + self.viewer_widget.mouse_origin
        )

        # Check for click selection (small displacement)
        dist_sq = self.current_dx * self.current_dx + self.current_dy * self.current_dy
        if dist_sq < 100.0:
            if hasattr(self.viewer_widget, "on_canvas_click"):
                self.viewer_widget.on_canvas_click(shift=shift)

        # Reset temporary drag delta
        self.viewer_widget.mouse_rotation = np.eye(3, dtype=float)
        self.viewer_widget.mouse_origin = np.zeros(3, dtype=float)
        self.is_dragging = False

        self.viewer_widget.request_render()

    def on_wheel(self, delta_y: float) -> None:
        """Handle mouse wheel zooming.

        Parameters
        ----------
        delta_y : float
            Wheel scroll delta (positive is scroll down/out, negative is scroll up/in).
        """
        # Browsers: deltaY > 0 indicates scrolling down/away (zoom out)
        # Using 2 ** (sensitivity * delta_y)
        zoom_factor = 2.0 ** (self.mouse_zoom_sensitivity * delta_y)
        self.viewer_widget.view_radius *= zoom_factor

        if self.viewer_widget.view_radius < self.min_view_radius:
            self.viewer_widget.view_radius = self.min_view_radius

        self.viewer_widget.request_render()

