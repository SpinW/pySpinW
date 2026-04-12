""" Representation of sites (e.g. magnetic atoms) within a magnetic system"""

import numpy as np
from numpy._typing import ArrayLike

from pyspinw.constants import ELECTRON_G
from pyspinw.serialisation import SPWSerialisationContext, SPWSerialisable, SPWDeserialisationContext, \
    numpy_deserialise, numpy_serialise

_unique_id_counter = -1
def _generate_unique_id():
    """ Generate a unique ID for each site currently loaded"""
    global _unique_id_counter # noqa: PLW0603
    _unique_id_counter += 1
    return _unique_id_counter

class LatticeSite(SPWSerialisable):
    """A spin site within a lattice

    :param: i,j,k - *required* Fractional coordinates within unit cell
    :param: sx, sy, sz - Spin components in the Cartesian directions defined by x||a, z perpendicular to a-b
    :param: supercell_spins - Spin for each propagation vector
    :param: g - g-tensor (3x3)
    :param: name
    """

    serialisation_name = "site"

    def __init__(self,
                 i: float, j: float, k: float,
                 sx: float | None = None,
                 sy: float | None = None,
                 sz: float | None = None,
                 supercell_spins: ArrayLike | None = None,
                 g: ArrayLike | None = None,
                 name: str = ""):

        self._i = float(i)
        self._j = float(j)
        self._k = float(k)

        #
        # Lots of case checking for the spin input format
        #

        if supercell_spins is None:
            self._spin_data = np.array([[
                0.0 if sx is None else sx,
                0.0 if sy is None else sy,
                0.0 if sz is None else sz]],
                    dtype=float)

        else:
            if sx is not None or sy is not None or sz is not None:
                raise ValueError("You need to specify either 'sx', 'sy' and 'sz', or 'supercell_spins'")

            good_shape = True

            supercell_spins = np.array(supercell_spins)

            if len(supercell_spins.shape) == 1:
                if supercell_spins.shape[0] != 3:
                    good_shape = False

            elif len(supercell_spins.shape) == 2:
                if supercell_spins.shape[1] != 3:
                    good_shape = False

            else:
                good_shape = False

            if not good_shape:
                raise ValueError("'supercell_spins' should have shape (3,) or (n,3)")

            self._spin_data = supercell_spins.reshape((-1, 3))

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

        self._base_spin = np.sum(self._spin_data, axis=0)
        self._name = name

        self._ijk = np.array([i, j, k], dtype=float)
        self._values = np.concatenate((self._ijk, self._base_spin))
        self._unique_id = _generate_unique_id()

    def n_components(self):
        """ Number of spin components for this site """
        return self.spin_data.shape[0]

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
    def base_spin(self):
        """spins as numpy array"""
        return self._base_spin

    @property
    def spin_data(self):
        """ Get all the spin data"""
        return self._spin_data

    @spin_data.setter
    def spin_data(self, spin_data):
        spin_data = np.array(spin_data)

        try:
            spin_data = spin_data.reshape(3, -1)

        except ValueError as e:
            raise ValueError("Expected spin data to be length 3, or convertable to a 3-by-n array")

        self._spin_data = spin_data % 1


    @property
    def g(self) -> np.ndarray:
        """ Magnetic g-factor, a 3x3 matrix/tensor """
        return self._g

    @property
    def values(self):
        """ ijk and spins as a numpy 6-vector"""
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
            sx=float(coordinates[3]),
            sy=float(coordinates[4]),
            sz=float(coordinates[5]),
            name=name)

    def __hash__(self):
        return self._unique_id

    def __repr__(self):
        m = self.base_spin

        if self.name is None or self.name == "":
            if np.sum(m**2) < 1e-9:
                return f"Site({self.i:.4g}, {self.j:.4g}, {self.k:.4g})"

            else:
                return f"Site({self.i:.4g}, {self.j:.4g}, {self.k:.4g}, spin={self.base_spin})"
        else:
            if np.sum(m ** 2) < 1e-9:
                return f"Site({self.name}, {self.i:.4g}, {self.j:.4g}, {self.k:.4g})"

            else:
                return f"Site({self.name}, {self.i:.4g}, {self.j:.4g}, {self.k:.4g}, spin={self.base_spin})"

    def _serialise(self, context: SPWSerialisationContext) -> dict:
        if not context.sites.has(self._unique_id):
            json = {
                "i": self.i,
                "j": self.j,
                "k": self.k,
                "supercell_spins": numpy_serialise(self._spin_data),
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
                    supercell_spins=numpy_deserialise(json["supercell_spins"]),
                    name=json["name"])
            else:
                site = LatticeSite(
                    i=json["i"],
                    j=json["j"],
                    k=json["k"],
                    supercell_spins=numpy_deserialise(json["supercell_spins"]),
                    name=json["name"])

            context.sites.put(response.id, site)

            return site


class ImpliedLatticeSite(LatticeSite):
    """ Lattice site that is implied by symmetry by a specified site"""

    def __init__(self,
                 parent_site: LatticeSite,
                 i: float, j: float, k: float,
                 sx: float | None = None,
                 sy: float | None = None,
                 sz: float | None = None,
                 supercell_spins: np.ndarray | None = None,
                 name: str | None = None):

        self._parent_site = parent_site

        super().__init__(i=i, j=j, k=k,
                         sx=sx, sy=sy, sz=sz,
                         supercell_spins=supercell_spins,
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
                "supercell_spins": numpy_serialise(self._spin_data),
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
            sx=float(coordinates[3]),
            sy=float(coordinates[4]),
            sz=float(coordinates[5]),
            name=name)

