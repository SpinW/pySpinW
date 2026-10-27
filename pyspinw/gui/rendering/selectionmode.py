"""Selection state of selectable objects."""

from enum import Enum


class SelectionMode(Enum):
    """Selection state of selectable objects.

    Attributes
    ----------
    NOT_SELECTED : int
        Object is neither selected nor hovered (``0``).
    HOVER : int
        Object is currently hovered over by the cursor (``1``).
    SELECTED : int
        Object is currently selected (``2``).
    SELECTED_HOVER : int
        Object is both selected and hovered over by the cursor (``3``).
    """

    NOT_SELECTED = 0
    HOVER = 1
    SELECTED = 2
    SELECTED_HOVER = 3

