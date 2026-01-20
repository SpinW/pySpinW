""" Rendering base classes """

import numpy as np

from pyspinw.gui.rendering.model import Model


class Renderer:
    """ Base class for rendering, renders objects of type pyspinw.gui.rendering.Model """

    def modelMatrix(self, rotation: np.ndarray, translation: np.ndarray, scale: np.ndarray):
        """ Get the model matrix in the GL MVP formalism from more intuitive values"""

    def render(self, model: Model, rotation=None, scaling=(1,1,1)):
        """ Render a model using this renderer"""

