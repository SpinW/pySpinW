import numpy as np
from pydantic import BaseModel


class CellOffset(BaseModel):
    """ Representation of the relative position of individual unit cells within the lattice """

    i: int
    j: int
    k: int

    _vector: np.ndarray | None = None

    def model_post_init(self, __context):
        self._vector = np.array([self.i, self.j, self.k], dtype=int)

    @property
    def as_tuple(self):
        return self.i, self.j, self.k

    @property
    def vector(self):
        return self._vector

    def position_in_supercell(self, supercell_size: tuple[int, int, int]):
        si, sj, sk = supercell_size

        return CellOffset(
            i=self.i % si,
            j=self.j % sj,
            k=self.k % sk)
