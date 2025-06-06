from typing import Annotated

import numpy as np

from pydantic import BaseModel, Field

_id_counter = -1
def _generate_unique_id():
    """ Generate a unique ID for each site - there must be a better way of doing this"""
    global _id_counter
    _id_counter += 1
    return _id_counter

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
    _unique_id: int | None = None

    @staticmethod
    def create(i: float, j: float, k: float,
               mi: float = 0, mj: float = 0, mk: float = 0,
               name: str | None = None):

        return LatticeSite(i=i, j=j, k=k, mi=mi, mj=mj, mk=mk, name=name)

    def model_post_init(self, __context):
        self._ijk = np.array([self.i, self.j, self.k], dtype=float)
        self._m = np.array([self.mi, self.mj, self.mk], dtype=float)
        self._values = np.concatenate((self._ijk, self._m))
        self._unique_id = _generate_unique_id()

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

    def __hash__(self):
        return self._unique_id

class ImpliedLatticeSite(LatticeSite):

    parent_site: LatticeSite

    @staticmethod
    def create(parent_site: LatticeSite,
               i: float, j: float, k: float,
               mi: float = 0, mj: float = 0, mk: float = 0,
               name: str | None = None):

        return ImpliedLatticeSite(parent_site=parent_site, i=i, j=j, k=k, mi=mi, mj=mj, mk=mk, name=name)

    @staticmethod
    def from_coordinates(parent_site: LatticeSite, coordinates: np.ndarray, name: str = ""):
        return ImpliedLatticeSite(
            parent_site=parent_site,
            i=coordinates[0],
            j=coordinates[1],
            k=coordinates[2],
            mi=coordinates[3],
            mj=coordinates[4],
            mk=coordinates[5],
            name=name)
