"""Different kinds of anisotropy"""
import numpy as np

from pyspinw.checks import check_sizes
from pyspinw.serialisation import SPWSerialisationContext, SPWSerialisable, numpy_serialise, numpy_deserialise, \
    SPWDeserialisationContext
from pyspinw.site import LatticeSite
from pyspinw.tolerances import tolerances


class Anisotropy(SPWSerialisable):
    """Defines the anisotropy at a given site"""

    @check_sizes(anisotropy_matrix=(3,3))
    def __init__(self, site: LatticeSite, anisotropy_matrix: np.ndarray):
        self._site = site
        self._anisotropy_matrix = anisotropy_matrix

    @property
    def anisotropy_matrix(self) -> np.ndarray:
        """Matrix spefifying the anisotropy - `A` term in the Hamiltonian"""
        if self._anisotropy_matrix is None:
            raise ValueError("Anisotropy matrix not initialised - this shouldn't happen")
        else:
            return self._anisotropy_matrix

    def _serialise(self, context: SPWSerialisationContext) -> dict:
        return {
            "site": self._site._serialise(context),
            "anisotropy_matrix": numpy_serialise(self._anisotropy_matrix)}

    @staticmethod
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        site = LatticeSite._deserialise(json["site"], context)
        anisotropy_matrix = numpy_deserialise(json["anisotropy_matrix"])
        return Anisotropy(site, anisotropy_matrix)


class DiagonalAnisotropy(Anisotropy):
    """Anisotropy oriented with axes, but variable amount in x, y and z"""

    @check_sizes(v=(3,), a=(1,), force_numpy=True)
    def __init__(self, site: LatticeSite, a: float, v: np.ndarray = np.array([0, 0, 1])):
        mag = np.sqrt(np.sum(v))

        if np.isclose(mag, 0, atol=tolerances.VECTOR_TOL):
            self._vector = np.zeros((3,), dtype=float)
        else:
            self._vector = v / mag

        anisotropy_matrix = a * v.reshape(3, 1) * v.reshape(1, 3)

        super().__init__(site, anisotropy_matrix)
        self._a = a

    @property
    def constant(self):
        return self._a

    @property
    def direction(self):
        return self._vector

    def _serialise(self, context: SPWSerialisationContext):
        return {
            "site": self._site._serialise(context),
            "vector": numpy_serialise(self._vector),
            "a": float(self._a)
        }

    @staticmethod
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        return DiagonalAnisotropy(
            LatticeSite._deserialise(json["site"], context),
            json["a"],
            numpy_deserialise(json["vector"]),
        )
