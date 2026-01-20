from numpy._typing import ArrayLike

from pyspinw.gui.rendering.model import Model
from pyspinw.gui.rendering.rendering import Renderable


class SelectionMode:
    NOT_SELECTED = 0
    SELECTED = 1
    SECONDARY = 2

class Site(Renderable):
    def __init__(self,
                 atom_model: Model,
                 moment_model: Model,
                 position: ArrayLike | None,
                 rotation: ArrayLike | None,
                 color: ArrayLike | None,
                 selection_mode: SelectionMode | int = SelectionMode.NOT_SELECTED,
                 moment_length: float = 1.5,
                 scale: float = 1.0,
                 anisotropy):

        self._moment_length = moment_length
        self._scale = scale

        super.__init__()

    @property
    def moment_length(self):
        pass

    @moment_length.setter
    def moment_length(self, moment_length):
        pass

    def scale(self):
        pass