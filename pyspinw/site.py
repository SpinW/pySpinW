""" Representation of sites (e.g. magnetic atoms) within a magnetic system"""

import numpy as np

from pydantic import BaseModel, Field

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
                 mi: float = 0.0, mj: float = 0.0, mk: float = 0.0,
                 name: str = ""):

        self._i = i
        self._j = j
        self._k = k

        self._mi = mi
        self._mj = mj
        self._mk = mk

        self._name = name

        self._ijk = np.array([i, j, k], dtype=float)
        self._m = np.array([mi, mj, mk], dtype=float)
        self._values = np.concatenate((self._ijk, self._m))
        self._unique_id = _generate_unique_id()

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
    def mi(self):
        """ Magnetic moment along first unit cell axis """
        return self._mi

    @property
    def mj(self):
        """ Magnetic moment along second unit cell axis """
        return self._mj

    @property
    def mk(self):
        """ Magnetic moment along third unit cell axis """
        return self._mk

    @property
    def ijk(self):
        """ ijk values as a numpy array"""
        return self._ijk

    @property
    def m(self):
        """magnetic moment as numpy array"""
        return self._m

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
        pass

    @staticmethod
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        site_or_json = context.sites.request_by_json(json)

        if site_or_json.deserialised:
            return site_or_json.value



class ImpliedLatticeSite(LatticeSite):
    """ Lattice site that is implied by symmetry by a specified site"""

    parent_site: LatticeSite

    @staticmethod
    def create(parent_site: LatticeSite,
               i: float, j: float, k: float,
               mi: float = 0, mj: float = 0, mk: float = 0,
               name: str | None = None):
        """ Create using ordered arguments"""
        return ImpliedLatticeSite(parent_site=parent_site, i=i, j=j, k=k, mi=mi, mj=mj, mk=mk, name=name)

    @staticmethod
    def from_coordinates(parent_site: LatticeSite, coordinates: np.ndarray, name: str = ""):
        """ Create ImpliedLatticeSite from coordinates"""
        return ImpliedLatticeSite(
            parent_site=parent_site,
            i=coordinates[0],
            j=coordinates[1],
            k=coordinates[2],
            mi=coordinates[3],
            mj=coordinates[4],
            mk=coordinates[5],
            name=name)

