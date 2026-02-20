""" Rendering base classes """

from abc import ABC, abstractmethod

from numpy._typing import ArrayLike
import numpy as np

from pyspinw.gui.camera import Camera
from pyspinw.gui.rendering.model import Model
from pyspinw.site import LatticeSite


class Renderable(ABC):
    def __init__(self,
                 position: ArrayLike | None,
                 rotation: ArrayLike | None,
                 scaling: ArrayLike | None,
                 color: ArrayLike | None,
                 ):

        # Set these fields to default and update, rather than calculate repeatedly
        self._position = np.zeros((3, ))
        self._rotation = np.eye(3)
        self._scaling = np.ones((3, ))
        self._color = np.ones((3, ))
        self._model_matrix = np.zeros((4, 4))

        self.position = position
        self.rotation = rotation
        self.scaling = scaling
        self.color = color

        self.update_model_matrix()

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, position: ArrayLike | None):
        self._position = np.zeros((1, 3)) if position is None \
            else np.array(position, dtype=np.float32).reshape((1, 3))

        self.update_model_matrix()

    @property
    def rotation(self):
        return self._rotation

    @rotation.setter
    def rotation(self, rotation: ArrayLike | None):
        self._rotation = np.eye(3, dtype=np.float32) if rotation is None \
            else np.array(rotation, dtype=np.float32)

        self.update_model_matrix()

    @property
    def scaling(self):
        return self._scaling

    @scaling.setter
    def scaling(self, scaling: ArrayLike | None):
        self._scaling = np.ones((3, ), dtype=np.float32) if scaling is None \
            else np.array(scaling, dtype=np.float32)

        self.update_model_matrix()


    @property
    def model_matrix(self):
        """ Get the model matrix in the GL MVP formalism from more intuitive values"""
        return self._model_matrix

    def update_model_matrix(self):
        scaling_matrix = np.array(np.diag(self._scaling), dtype=np.float32)

        filled_part = np.concatenate((self._rotation @ scaling_matrix, self._position), axis=1)

        self._model_matrix = np.concatenate((filled_part, np.array([0,0,0,1], dtype=np.float32)), axis=0)

if __name__ == "__main__":
    pass
