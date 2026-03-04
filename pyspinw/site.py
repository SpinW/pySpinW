""" Representation of sites (e.g. magnetic atoms) within a magnetic system"""

import numpy as np
from numpy._typing import ArrayLike
from scipy.stats import goodness_of_fit

from pyspinw.constants import ELECTRON_G
from pyspinw.serialisation import SPWSerialisationContext, SPWSerialisable, SPWDeserialisationContext, \
    numpy_deserialise, numpy_serialise

_id_counter = -1
def _generate_unique_id():
    """ Generate a unique ID for each site currently loaded"""
    global _id_counter # noqa: PLW0603
    _id_counter += 1
    return _id_counter

class LatticeSite(SPWSerialisable):
    """A spin site within a lattice

    :param: i,j,k - Fractional coordinates within unit cell
    :param: mi,mj,mk - Magnetic moment along unit cell aligned axis
    """

    serialisation_name = "site"

    def __init__(self,
                 i: float, j: float, k: float,
                 mi: float | None = None,
                 mj: float | None = None,
                 mk: float | None = None,
                 supercell_moments: ArrayLike | None = None,
                 g: ArrayLike | None = None,
                 name: str = ""):

        self._i = float(i)
        self._j = float(j)
        self._k = float(k)

        #
        # Lots of case checking for the moment input format
        #

        if supercell_moments is None:
            self._moment_data = np.array([[
                0.0 if mi is None else mi,
                0.0 if mj is None else mj,
                0.0 if mk is None else mk]],
                    dtype=float)

        else:
            if mi is not None or mj is not None or mk is not None:
                raise ValueError("You need to specify either 'mi', 'mj' and 'mk', or 'supercell_moments'")

            good_shape = True

            supercell_moments = np.array(supercell_moments)

            if len(supercell_moments.shape) == 1:
                if supercell_moments.shape[0] != 3:
                    good_shape = False

            elif len(supercell_moments.shape) == 2:
                if supercell_moments.shape[1] != 3:
                    good_shape = False

            else:
                good_shape = False

            if not good_shape:
                raise ValueError("'supercell_moments' should have shape (3,) or (n,3)")

            self._moment_data = supercell_moments.reshape((-1, 3))

        # Get the g tensor
        if g is None:
            self._g = ELECTRON_G * np.eye(3)
        else:
            g = np.array(g)

            flattened = g.reshape(-1)
            if flattened.shape[0] == 1:
                self._g = g * np.eye(3)
            elif flattened.shape[0] == 3:
                self._g = np.diag(flattened)
            elif flattened.shape[0] == 9:
                self._g = g
            else:
                raise ValueError("g-factor should be a scalar, a vector of length 3, or a 3-by-3 matrix")

        self._base_moment = np.sum(self._moment_data, axis=0)

        self._name = name

        self._ijk = np.array([i, j, k], dtype=float)
        self._values = np.concatenate((self._ijk, self._base_moment))
        self._unique_id = _generate_unique_id()

    @property
    def unique_id(self):
        """ Unique ID for this site """
        return self._unique_id

    @property
    def name(self):
        """ Name given to this site """
        return self._name

    @property
    def i(self):
        """ Fractional position along first unit cell axis """
        return self._i

    @property
    def j(self):
        """ Fractional position along second unit cell axis """
        return self._j

    @property
    def k(self):
        """ Fractional position along third unit cell axis """
        return self._k

    @property
    def ijk(self):
        """ ijk values as a numpy array"""
        return self._ijk

    @property
    def base_moment(self):
        """magnetic moment as numpy array"""
        return self._base_moment

    @property
    def moment_data(self):
        """ Get all the magnetic moment data"""
        return self._moment_data

    @property
    def g(self) -> np.ndarray:
        """ Magnetic g-factor, a 3x3 matrix/tensor """
        return self._g

    @property
    def values(self):
        """ ijk and moments as a numpy 6-vector"""
        return self._values

    @property
    def parent_site(self):
        """ Get the parent site (just itself for non-implied sites)"""
        return self

    @staticmethod
    def from_coordinates(coordinates: np.ndarray, name: str = ""):
        """ Create from an array of values """
        return LatticeSite(
            i=float(coordinates[0]),
            j=float(coordinates[1]),
            k=float(coordinates[2]),
            mi=float(coordinates[3]),
            mj=float(coordinates[4]),
            mk=float(coordinates[5]),
            name=name)

    def __hash__(self):
        return self._unique_id

    def __repr__(self):
        m = self.base_moment
        if np.sum(m**2) < 1e-9:
            return f"Site({self.i:.4g}, {self.j:.4g}, {self.k:.4g})"

        else:
            return f"Site({self.i:.4g}, {self.j:.4g}, {self.k:.4g}, base_moment={self.base_moment})"

    def _serialise(self, context: SPWSerialisationContext) -> dict:
        if not context.sites.has(self._unique_id):
            json = {
                "i": self.i,
                "j": self.j,
                "k": self.k,
                "supercell_moments": numpy_serialise(self._moment_data),
                "name": self.name,
                "g": numpy_serialise(self.g)
            }

            context.sites.put(self._unique_id, json)

        return context.sites.reference(self._unique_id)


    @staticmethod
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        response = context.sites.request_by_json(json)

        if response.deserialised:
            return response.value

        else:
            json = response.value
            # Do the actual deserialisation

            if "parent" in json:
                parent = LatticeSite._deserialise(json["parent"], context)
                site = ImpliedLatticeSite(
                    parent_site = parent,
                    i=json["i"],
                    j=json["j"],
                    k=json["k"],
                    supercell_moments=numpy_deserialise(json["supercell_moments"]),
                    name=json["name"])
            else:
                site = LatticeSite(
                    i=json["i"],
                    j=json["j"],
                    k=json["k"],
                    supercell_moments=numpy_deserialise(json["supercell_moments"]),
                    name=json["name"])

            context.sites.put(response.id, site)

            return site


class ImpliedLatticeSite(LatticeSite):
    """ Lattice site that is implied by symmetry by a specified site"""

    def __init__(self,
                 parent_site: LatticeSite,
                 i: float, j: float, k: float,
                 mi: float | None = None,
                 mj: float | None = None,
                 mk: float | None = None,
                 supercell_moments: np.ndarray | None = None,
                 name: str | None = None):

        self._parent_site = parent_site

        super().__init__(i=i, j=j, k=k,
                         mi=mi, mj=mj, mk=mk,
                         supercell_moments=supercell_moments,
                         name=name)

    @property
    def parent_site(self):
        """ Get the parent of this site """
        return self._parent_site

    def _serialise(self, context: SPWSerialisationContext) -> dict:
        if not context.sites.has(self._unique_id):
            parent_ref = self.parent_site._serialise(context)
            json = {
                "i": self.i,
                "j": self.j,
                "k": self.k,
                "supercell_moments": numpy_serialise(self._moment_data),
                "name": self.name,
                "parent": parent_ref
            }

            context.sites.put(self._unique_id, json)

        return context.sites.reference(self._unique_id)



    @staticmethod
    def from_coordinates(parent_site: LatticeSite, coordinates: np.ndarray, name: str = ""):
        """ Create ImpliedLatticeSite from coordinates"""
        return ImpliedLatticeSite(
            parent_site=parent_site,
            i=float(coordinates[0]),
            j=float(coordinates[1]),
            k=float(coordinates[2]),
            mi=float(coordinates[3]),
            mj=float(coordinates[4]),
            mk=float(coordinates[5]),
            name=name)

