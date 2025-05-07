"""Different kinds of anisotropy"""
import numpy as np

from pyspinw._base import Anisotropy, Identifier
from pyspinw.checks import check_sizes


class GeneralAnisotropy(Anisotropy):
    """General anisotropy specification"""

    @check_sizes(a=(3,3), force_numpy=True)
    def __init__(self, site: Identifier, a: np.ndarray):
        super().__init__(site)

        self._a = a

    @property
    def anisotropy_matrix(self) -> np.ndarray:
        """The matrix defining the anisotropy."""
        return self._a


class DiagonalAnisotropy(GeneralAnisotropy):
    """Anisotropy oriented with axes, but variable amount in x, y and z"""

    @check_sizes(a=(3,), force_numpy=True)
    def __init__(self, site: Identifier, a: np.ndarray):
        super().__init__(site, a)
        self._v = a


class XAxisAnisotropy(DiagonalAnisotropy):
    """Pure X anisotropy"""

    def __init__(self, site: Identifier, a: float):
        super().__init__(site, np.array([1,0,0], dtype=float))
        self.a_x = a


class YAxisAnisotropy(DiagonalAnisotropy):
    """Pure Y anisotropy"""

    def __init__(self, site: Identifier, a: float):
        super().__init__(site, np.array([0,1,0], dtype=float))
        self.a_y = a


class ZAxisAnisotropy(DiagonalAnisotropy):
    """Pure Z anisotropy"""

    def __init__(self, site: Identifier, a: float):
        super().__init__(site, np.array([0,0,1], dtype=float))
        self.a_z = a
