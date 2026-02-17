from numpy._typing import ArrayLike

from pyspinw.gui.camera import Camera
from pyspinw.gui.rendering.model import Model
from pyspinw.gui.rendering.rendering import Renderable


class SelectionMode:
    NOT_SELECTED = 0
    HOVER = 1
    SELECTED = 2
    SELECTED_HOVER = 3

