import numpy as np

from pyspinw.gui.rendering.model import Model


class Renderer:
    def modelMatrix(self, rotation: np.ndarray, translation: np.ndarray, scale: np.ndarray):
        pass


    def render(self, model: Model, rotation=None, scaling=(1,1,1)):
        pass

