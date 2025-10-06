""" Representation of sites (e.g. magnetic atoms) within a magnetic system"""

import numpy as np
from scipy.stats import goodness_of_fit

from pyspinw.serialisation import SPWSerialisationContext, SPWSerialisable, SPWDeserialisationContext

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
                 supercell_moments: np.ndarray | None = None,
                 name: str = ""):

        self._i = i
        self._j = j
        self._k = k

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
                raise ValueError("You need to specify at least one of 'mi', 'mj', 'mk' or 'supercell_moments'")

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

        self._base_moment = np.sum(supercell_moments, axis=0)

        self._name = name

        self._ijk = np.array([i, j, k], dtype=float)
        self._m = np.array([mi, mj, mk], dtype=float)
        self._values = np.concatenate((self._ijk, self._m))
        self._unique_id = _generate_unique_id()

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
        return self._m

    @property
    def moment_data(self):
        """ Get all the magnetic moment data"""
        return self._moment_data

    @property
    def values(self):
        """ ijk and moments as a numpy 6-vector"""
        return self._values

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

    def _serialise(self, context: SPWSerialisationContext) -> dict:
        if not context.sites.has(self._unique_id):
            json = {
                "i": self.i,
                "j": self.j,
                "k": self.k,
                "mi": self.mi,
                "mj": self.mj,
                "mk": self.mk,
                "name": self.name
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
                    mi=json["mi"],
                    mj=json["mj"],
                    mk=json["mk"],
                    name=json["name"])
            else:
                site = LatticeSite(
                    i=json["i"],
                    j=json["j"],
                    k=json["k"],
                    mi=json["mi"],
                    mj=json["mj"],
                    mk=json["mk"],
                    name=json["name"])

            context.sites.put(response.id, site)

            return site


class ImpliedLatticeSite(LatticeSite):
    """ Lattice site that is implied by symmetry by a specified site"""

    def __init__(self,
                 parent_site: LatticeSite,
                 i: float, j: float, k: float,
                 mi: float = 0, mj: float = 0, mk: float = 0,
                 name: str | None = None):

        self.parent_site = parent_site

        super().__init__(i=i, j=j, k=k, mi=mi, mj=mj, mk=mk, name=name)


    def _serialise(self, context: SPWSerialisationContext) -> dict:
        if not context.sites.has(self._unique_id):
            parent_ref = self.parent_site._serialise(context)
            json = {
                "i": self.i,
                "j": self.j,
                "k": self.k,
                "mi": self.mi,
                "mj": self.mj,
                "mk": self.mk,
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

