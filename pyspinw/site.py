from typing import Annotated

import numpy as np

from pydantic import BaseModel, Field

from pyspinw.symmetry.unitcell import UnitCell


class LatticeSite(BaseModel):
    """A spin site within a lattice """

    i: Annotated[float, Field(ge=0.0, lt=1.0)]
    j: Annotated[float, Field(ge=0.0, lt=1.0)]
    k: Annotated[float, Field(ge=0.0, lt=1.0)]

    mi: float = 0.0
    mj: float = 0.0
    mk: float = 0.0

    ijk: float | None = Field(default=None, init=False)
    m: float | None = Field(default=None, init=False)


    def __init__(self,
                 i: float, j: float, k: float,
                 mi: float = 0, mj: float = 0, mk: float = 0):

        super().__init__(i=i, j=j, k=k, mi=mi, m=mj, mk=mk)

    def model_post_init(self, __context):
        self.ijk = np.array([self.i, self.j, self.k], dtype=float)
        self.m = np.array([self.mi, self.mj, self.mk], dtype=float)
