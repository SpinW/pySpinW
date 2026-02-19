"""Different kinds of anisotropy"""
import numpy as np
from numpy._typing import ArrayLike

from pyspinw.checks import check_sizes
from pyspinw.serialisation import SPWSerialisationContext, SPWSerialisable, numpy_serialise, numpy_deserialise, \
    SPWDeserialisationContext
from pyspinw.site import LatticeSite
from pyspinw.tolerances import tolerances


class Anisotropy(SPWSerialisable):
    """Defines the anisotropy at a given site"""

    serialisation_name = "anisotropy"

    @check_sizes(anisotropy_matrix=(3,3))
    def __init__(self, site: LatticeSite, anisotropy_matrix: ArrayLike):
        self._site = site
        self._anisotropy_matrix = np.array(anisotropy_matrix)

    @property
    def site(self):
        """ Get the site for this anisotropy """
        return self._site

    @property
    def anisotropy_matrix(self) -> np.ndarray:
        """Matrix specifying the anisotropy - `A` term in the Hamiltonian"""
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

    @property
    def parameter_string(self):
        m = self.anisotropy_matrix.reshape(-1)
        return f"[[{m[0]}, {m[1]}, {m[2]}], [{m[3]}, {m[4]}, {m[5]}], [{m[6]}, {m[7]}, {m[8]}]]"

    def __repr__(self):
        return (f"Anisotropy({self.site.name}, {self.parameter_string})")

    def updated(self, site: LatticeSite | None = None, anisotropy_matrix: ArrayLike | None = None):
        """ Return a copy of this anisotropy term with variables replaced"""
        return Anisotropy(
            site=self.site if site is None else site,
            anisotropy_matrix = self.anisotropy_matrix if anisotropy_matrix is None else np.array(anisotropy_matrix))

class AxisMagnitudeAnisotropy(Anisotropy):
    """Anisotropy oriented with axes, but variable amount in x, y and z"""

    @check_sizes(direction=(3,), force_numpy=True)
    def __init__(self, site: LatticeSite, a: float, direction: np.ndarray = np.array([0, 0, 1])):
        mag = np.sqrt(np.sum(direction**2))

        if np.isclose(mag, 0, atol=tolerances.VECTOR_TOL):
            self._direction = np.zeros((3,), dtype=float)
        else:
            self._direction = direction / mag

        anisotropy_matrix = a * direction.reshape(3, 1) * direction.reshape(1, 3)

        super().__init__(site, anisotropy_matrix)
        self._a = a

    @property
    def a(self):
        """ Amount of anisotropy (anisotropy constant)"""
        return self._a

    @property
    def constant(self):
        """ Amount of anisotropy (anisotropy constant), alias for `a` """
        return self._a

    @property
    def direction(self):
        """ The principal direction of the anisotropy"""
        return self._direction

    @property
    def parameter_string(self):
        return f"a={self.constant}, axis={self.direction}"

    def __repr__(self):
        return (f"Anisotropy({self.site.name}, {self.parameter_string})")

    def _serialise(self, context: SPWSerialisationContext):
        return {
            "site": self._site._serialise(context),
            "direction": numpy_serialise(self._direction),
            "a": float(self._a)
        }

    @staticmethod
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        return AxisMagnitudeAnisotropy(
            LatticeSite._deserialise(json["site"], context),
            json["a"],
            numpy_deserialise(json["direction"]),
        )

    def updated(self, site: LatticeSite | None = None, a: float | None = None, direction: ArrayLike | None = None):
        """ Return a copy of this anisotropy term with variables replaced"""
        return AxisMagnitudeAnisotropy(
            site=self.site if site is None else site,
            a=self.a if a is None else a,
            direction=self.direction if direction is None else np.array(direction))
