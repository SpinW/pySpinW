""" Cell offsets (i.e. reference to a particular unit cell in a crystal)"""

import numpy as np
from pydantic import BaseModel


class CellOffset(BaseModel):
    """ Representation of the relative position of individual unit cells within the lattice """

    i: int
    j: int
    k: int

    _vector: np.ndarray | None = None

    def model_post_init(self, __context):
        """pydantic: after init"""
        self._vector = np.array([self.i, self.j, self.k], dtype=int)

    @property
    def as_tuple(self):
        """ As a tuple of ints"""
        return self.i, self.j, self.k

    @property
    def vector(self):
        """ Vector of index (int array)"""
        return self._vector

    def position_in_supercell(self, supercell_size: tuple[int, int, int]):
        """ Get the position within a supercell, rather than an infinite crystal (do mods)"""
        si, sj, sk = supercell_size

        return CellOffset(
            i=self.i % si,
            j=self.j % sj,
            k=self.k % sk)
