""" Supercell representations"""
from abc import ABC, abstractmethod
from fractions import Fraction
from math import lcm

import numpy as np

from pyspinw.serialisation import SPWSerialisable, SPWSerialisationContext, \
    serialise_fraction_or_builtin, deserialise_fraction_or_builtin, SPWDeserialisationContext, expects_keys
from pyspinw.site import LatticeSite
from pyspinw.cell_offsets import CellOffset
from pyspinw.util import rotation_matrix


def _coerce_numeric_input(value: Fraction | float, max_denom=1000_000):
    """ Convert floats or fractions into the form needed for a propagation vector"""
    if not isinstance(value, Fraction):
        value = Fraction(value).limit_denominator(max_denominator=max_denom)

    # It's conventional to use 0 rather than one for a period 1 supercell (as they're indistinguishable)
    # but mathematically, it should be one

    if value == 0:
        return Fraction(1)

    else:
        return value


class PropagationVector(SPWSerialisable):
    """ Propagation vector"""

    serialisation_name = "propagation_vector"

    def __init__(self,
                i: Fraction | float,
                j: Fraction | float,
                k: Fraction | float):

        self.i = i
        self.j = j
        self.k = k

        self._vector = np.array([self.i, self.j, self.k], dtype=float)


    @property
    def vector(self):
        """ As a numpy vector"""
        return self._vector

    def dot(self, cell_offset: CellOffset):
        """ Dot product with a cell offset"""
        return np.dot(self.vector, cell_offset.vector)

    def _serialise(self, context: SPWSerialisationContext) -> dict:
        return {
            "i": serialise_fraction_or_builtin(self.i),
            "j": serialise_fraction_or_builtin(self.j),
            "k": serialise_fraction_or_builtin(self.k),
            "assured_commensurate": False
        }

    @staticmethod
    @expects_keys("i,j,k")
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        i = deserialise_fraction_or_builtin(json["i"])
        j = deserialise_fraction_or_builtin(json["j"])
        k = deserialise_fraction_or_builtin(json["k"])

        if json["assured_commensurate"]:
            return CommensuratePropagationVector(i, j, k)
        else:
            return PropagationVector(i, j, k)


class CommensuratePropagationVector(PropagationVector):
    """ Propagation vector with integer values"""

    def __init__(self,
                 i: Fraction | float,
                 j: Fraction | float,
                 k: Fraction | float):

        i = _coerce_numeric_input(i)
        j = _coerce_numeric_input(j)
        k = _coerce_numeric_input(k)

        super().__init__(i, j, k)

        self._vector = np.array([self.i, self.j, self.k], dtype=float)


    def _serialise(self, context: SPWSerialisationContext) -> dict:
        serialised = super()._serialise(context)
        serialised["assured_commensurate"] = True
        return serialised

class SupercellTransformation(ABC, SPWSerialisable):

    serialisation_name = "supercell_transformation"

    @abstractmethod
    def apply(self, moment: np.ndarray, propagation_vector: CommensuratePropagationVector, cell_offset: CellOffset):
        """ Apply this transformation to a moment with a given offset """



class IdentityTransform(SupercellTransformation):
    def apply(self, moment: np.ndarray, propagation_vector: CommensuratePropagationVector, cell_offset: CellOffset):
        return moment

class SupercellRotation(SupercellTransformation):
    """ Rotation transformation in a supercell"""
    def __init__(self, axis: np.ndarray):
        self._axis = axis

    def apply(self, moment: np.ndarray, propagation_vector: CommensuratePropagationVector, cell_offset: CellOffset):
        angle = 2 * np.pi * propagation_vector.dot(cell_offset)
        return moment @ rotation_matrix(angle, self._axis)


class Supercell(ABC, SPWSerialisable):

    serialisation_name = "supercell"

    def __init__(self, scaling: tuple[int, int, int] | None = None):

        if any([scaling_axis < 1 for scaling_axis in scaling]):
            raise ValueError("Scaling components should be >= 1")

        if scaling is None:
            self._scaling = (1,1,1)
        else:
            self._scaling = scaling

    @abstractmethod
    def moment(self, site: LatticeSite, cell_offset: CellOffset):
        """ Evaluate the moment of a lattice site in this unit cell """

    @abstractmethod
    def cell_size(self) -> tuple[int, int, int]:
        """ How big is this supercell """

    @abstractmethod
    def summation_form(self) -> "Supercell":
        """Get a summation type supercell"""

    def _serialise(self, context: SPWSerialisationContext) -> dict:
        pass

    def _deserialise(json: dict, context: SPWDeserialisationContext):
        pass


class NoSupercell(Supercell):
    """ Trivial supercell, just a single unit cell """
    def moment(self, site: LatticeSite, cell_offset: CellOffset):
        return site.base_moment

    def cell_size(self) -> tuple[int, int, int]:
        return self._scaling

    def summation_form(self) -> "Supercell":
        return self


class CommensurateSupercell(Supercell):
    """ Base class for a commensurate supercell"""

    def __init__(self, propagation_vectors, scaling: tuple[int, int, int] | None = None):
        self._propagation_vectors = propagation_vectors
        super().__init__(scaling)


    def _transform_evaluate(self, cell_offset: CellOffset, moment: np.ndarray):
        """ Apply the transformation to a given moment """
        raise NotImplementedError("evaluate not implemented in base class")

    def moment(self, site: LatticeSite, cell_offset: CellOffset):
        return self._transform_evaluate(cell_offset=cell_offset, moment=site.base_moment)

    def cell_size(self) -> tuple[int, int, int] | None:
        """ Get the smallest possible supercell"""
        i = lcm(*[vector.i.denominator for vector in self.propagation_vectors])
        j = lcm(*[vector.j.denominator for vector in self.propagation_vectors])
        k = lcm(*[vector.k.denominator for vector in self.propagation_vectors])

        return self._scaling[0]*i, self._scaling[1]*j, self._scaling[2]*k

    def n_cells(self):
        """ Number of cells in this supercell """
        a, b, c = self.cell_size()
        return a*b*c


class TransformationSupercell(CommensurateSupercell):
    """ Supercell with moments defined by the following equation:

    m = prod(M_i^{r.k}) m0

    where m0 is the existing moment at a given site
    """

    def __init__(self, *transforms: tuple[CommensuratePropagationVector, SupercellTransformation | None]):
        # TODO: Provide a nicer interface for this maybe
        self._transforms = [(vector, IdentityTransform() if transform is None else transform)
                            for vector, transform in transforms]


    def _transform_evaluate(self, cell_offset: CellOffset, moment: np.ndarray):
        for vector, transform in self._transforms:
            moment = transform.apply(moment=moment, propagation_vector=vector, cell_offset=cell_offset)

        return moment

    def summation_form(self) -> "Supercell":
        pass


class SummationSupercell(CommensurateSupercell):
    """ Supercell with moments defined by

    m = sum(m_j exp(2 pi i r.k_j)).real
    """

    def __init__(self, propagation_vectors):
        super().__init__(propagation_vectors)


    def moment(self, cell_offset: CellOffset, moment: np.ndarray):
        """ Calculate moment at a given cell offset"""
        return np.sum(component.evaluate(cell_offset) for component in self.components).real

    def summation_form(self, moment) -> "SummationSupercell":
        """ Convert to summation form (it is already in this form, but not all Supercells are)"""
        return self
