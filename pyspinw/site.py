from typing import Annotated

import numpy as np

from pydantic import BaseModel, Field


class LatticeSite(BaseModel):
    """A spin site within a lattice """

    i: float
    j: float
    k: float

    mi: float = 0.0
    mj: float = 0.0
    mk: float = 0.0

    name: str | None = None

    _ijk: np.ndarray | None = None
    _m: np.ndarray | None = None
    _values: np.ndarray | None = None

    def __init__(self,
                 i: float, j: float, k: float,
                 mi: float = 0, mj: float = 0, mk: float = 0, name: str | None = None):

        super().__init__(i=i, j=j, k=k, mi=mi, mj=mj, mk=mk, name=name)

    def model_post_init(self, __context):
        self._ijk = np.array([self.i, self.j, self.k], dtype=float)
        self._m = np.array([self.mi, self.mj, self.mk], dtype=float)
        self._values = np.concatenate((self._ijk, self._m))

    @property
    def ijk(self):
        return self._ijk

    @property
    def m(self):
        return self._m

    @property
    def values(self):
        return self._values

    @staticmethod
    def from_coordinates(coordinates: np.ndarray, name: str = ""):
        return LatticeSite(
            i=coordinates[0],
            j=coordinates[1],
            k=coordinates[2],
            mi=coordinates[3],
            mj=coordinates[4],
            mk=coordinates[5],
            name=name)