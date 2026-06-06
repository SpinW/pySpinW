"""Different kinds of anisotropy."""

import numpy as np
from numpy._typing import ArrayLike

from pyspinw.checks import check_sizes
from pyspinw.serialisation import (
    SPWSerialisationContext,
    SPWSerialisable,
    numpy_serialise,
    numpy_deserialise,
    SPWDeserialisationContext,
)
from pyspinw.site import LatticeSite
from pyspinw.tolerances import tolerances


_anisotropy_id_counter = -1


def _generate_unique_anisotropy_id():
    """Generate a unique ID for each anisotropy currently loaded."""
    global _anisotropy_id_counter  # noqa: PLW0603
    _anisotropy_id_counter += 1
    return _anisotropy_id_counter


class Anisotropy(SPWSerialisable):
    """Represent a single-site anisotropy term.

    Parameters
    ----------
    site
        Lattice site to which the anisotropy term applies.
    anisotropy_matrix
        3x3 matrix defining the anisotropy contribution for the site.
    """

    #: Type name used when serialising anisotropy terms to SPW data.
    serialisation_name = "anisotropy"

    #: Names of scalar fields that can be varied by parameter definitions.
    scalar_parameters = []

    @check_sizes(anisotropy_matrix=(3, 3))
    def __init__(self, site: LatticeSite, anisotropy_matrix: ArrayLike):
        self._site = site
        self._anisotropy_matrix = np.array(anisotropy_matrix)
        self._unique_id = _generate_unique_anisotropy_id()

    @property
    def unique_id(self):
        """Unique ID for this anisotropy."""
        return self._unique_id

    @property
    def site(self):
        """Get the site for this anisotropy."""
        return self._site

    @property
    def anisotropy_matrix(self) -> np.ndarray:
        """Matrix specifying the anisotropy - `A` term in the Hamiltonian."""
        if self._anisotropy_matrix is None:
            raise ValueError("Anisotropy matrix not initialised - this shouldn't happen")
        else:
            return self._anisotropy_matrix

    def _serialise(self, context: SPWSerialisationContext) -> dict:
        return {"site": self._site._serialise(context), "anisotropy_matrix": numpy_serialise(self._anisotropy_matrix)}

    @staticmethod
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        site = LatticeSite._deserialise(json["site"], context)
        anisotropy_matrix = numpy_deserialise(json["anisotropy_matrix"])
        return Anisotropy(site, anisotropy_matrix)

    @property
    def parameter_string(self):
        """A string representation of the parameters."""
        m = self.anisotropy_matrix.reshape(-1)
        return f"[[{m[0]}, {m[1]}, {m[2]}], [{m[3]}, {m[4]}, {m[5]}], [{m[6]}, {m[7]}, {m[8]}]]"

    def __repr__(self):
        return f"Anisotropy({self.site.name}, {self.parameter_string})"

    def updated(self, site: LatticeSite | None = None, anisotropy_matrix: ArrayLike | None = None):
        """Return a copy of this anisotropy term with variables replaced.

        Parameters
        ----------
        site
            Replacement lattice site. If omitted, the current site is reused.
        anisotropy_matrix
            Replacement anisotropy matrix. If omitted, the current matrix is reused.
        """
        return Anisotropy(
            site=self.site if site is None else site,
            anisotropy_matrix=self.anisotropy_matrix if anisotropy_matrix is None else np.array(anisotropy_matrix),
        )


class AxisMagnitudeAnisotropy(Anisotropy):
    """Represent a uniaxial anisotropy term.

    Parameters
    ----------
    site
        Lattice site to which the anisotropy term applies.
    a
        Anisotropy constant.
    direction
        Principal anisotropy direction. The vector is normalized before use.
    """

    #: Names of scalar fields that can be varied by parameter definitions.
    scalar_parameters = ["a"]

    @check_sizes(direction=(3,), force_numpy=True)
    def __init__(self, site: LatticeSite, a: float, direction: ArrayLike = np.array([0, 0, 1])):
        direction = np.array(direction)
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
        """Amount of anisotropy (anisotropy constant)."""
        return self._a

    @property
    def constant(self):
        """Amount of anisotropy (anisotropy constant), alias for `a`."""
        return self._a

    @property
    def direction(self):
        """The principal direction of the anisotropy."""
        return self._direction

    @property
    def parameter_string(self):
        """A string representation of the parameters."""
        return f"a={self.constant}, axis={self.direction}"

    def __repr__(self):
        return f"Anisotropy({self.site.name}, {self.parameter_string})"

    def _serialise(self, context: SPWSerialisationContext):
        return {
            "site": self._site._serialise(context),
            "direction": numpy_serialise(self._direction),
            "a": float(self._a),
        }

    @staticmethod
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        return AxisMagnitudeAnisotropy(
            LatticeSite._deserialise(json["site"], context),
            json["a"],
            numpy_deserialise(json["direction"]),
        )

    def updated(self, site: LatticeSite | None = None, a: float | None = None, direction: ArrayLike | None = None):
        """Return a copy of this anisotropy term with variables replaced.

        Parameters
        ----------
        site
            Replacement lattice site. If omitted, the current site is reused.
        a
            Replacement anisotropy constant. If omitted, the current constant is reused.
        direction
            Replacement principal anisotropy direction. If omitted, the current direction is reused.
        """
        return AxisMagnitudeAnisotropy(
            site=self.site if site is None else site,
            a=self.a if a is None else a,
            direction=self.direction if direction is None else np.array(direction),
        )
