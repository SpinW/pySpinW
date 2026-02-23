""" Selection state of selectable objects"""

from enum import Enum


class SelectionMode(Enum):
    """ Selection state of selectable objects """

    NOT_SELECTED = 0
    HOVER = 1
    SELECTED = 2
    SELECTED_HOVER = 3

