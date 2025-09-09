""" Supercell representations"""

from fractions import Fraction
from math import lcm
from typing import Any

import numpy as np

from checks import check_sizes
from pyspinw.tolerances import tolerances
from pyspinw.cell_offsets import CellOffset

class PropagationVector:
    """ Propagation vector"""

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

class IncommensuratePropagationVector(PropagationVector):
    """ Propagation vector with non-integral values"""


class CommensurateSupercell:
    """ Base class for a commensurate supercell"""

    def __init__(self, propagation_vectors):
        propagation_vectors = propagation_vectors

    def evaluate(self, cell_offset: CellOffset, moment: np.ndarray):
        """ Get the moment at a given cell location"""
        raise NotImplementedError("evaluate not implemented in base class")

    def minimal_supercell(self) -> tuple[int, int, int] | None:
        """ Get the smallest possible supercell"""
        i = lcm(*[vector.i.denominator for vector in self.propagation_vectors])
        j = lcm(*[vector.j.denominator for vector in self.propagation_vectors])
        k = lcm(*[vector.k.denominator for vector in self.propagation_vectors])

        return i, j, k

    def summation_form(self, moment) -> "SummationSupercell":
        """Get a summation type supercell"""
        raise NotImplementedError("summation_form not implemented in base class")

class SummationSupercell(CommensurateSupercell):
    """ Supercell with moments defined by

    m = sum(m_j exp(2 pi i r.k_j)).real
    """

    def __init__(self, propagation_vectors):
        super().__init__(propagation_vectors)


    def apply(self, moments):
        for moment in moments:
            if not isinstance(moment, np.ndarray):
                raise TypeError("moments must contain numpy ndarrays")
            if moment.shape != (3,):
                raise ValueError("moments must be of shape (3,)")

    def evaluate(self, cell_offset: CellOffset, moment: np.ndarray):
        """ Calculate moment at a given cell offset"""
        return np.sum(component.evaluate(cell_offset) for component in self.components).real

    def summation_form(self, moment) -> "SummationSupercell":
        """ Convert to summation form (it is already in this form, but not all Supercells are)"""
        return self

class SupercellTransformation:

    @check_sizes(matrix=(3,3))
    def __init__(self, matrix: np.ndarray):
        self._matrix = matrix

    @property
    def matrix(self):
        return self._matrix

class GeneralSupercellTransformation

class SupercellRotation(SupercellTransformation):
    def __init__(self, propagation_vector: CommensuratePropagationVector, axis: np.ndarray):




class TransformationSupercell(CommensurateSupercell):
    """ Supercell with moments defined by the following equation:

    m = prod(M_i^{r.k}) m0

    where m0 is the existing moment at a given site
    """

    transforms: list[np.ndarray]

    @model_validator(mode="after")
    def check_lengths(self):
        """pydantic: checks sizes"""
        if len(self.transforms) != len(self.propagation_vectors):
            raise ValueError('the lists: transforms and propagation_vectors must have the same length')
        return self



    def evaluate(self, cell_offset: CellOffset, moment: np.ndarray):
        """ Get the moment at a given cell offset"""
        # Get the values of k.r, reduce to unit cell though, which will make sure we have integer powers
        powers = [propagation_vector.dot(cell_offset.position_in_supercell(self.minimal_supercell()))
                  for propagation_vector in self.propagation_vectors]

        int_powers = [int(power) for power in powers]

        # Check they're really integers
        for power, int_power in zip(powers, int_powers):
            if abs(power - int_power) > tolerances.IS_INTEGER_TOL:
                raise ValueError(f"Supercell k.r has non-integer values ({power})")

        # TODO: Certainly better performing ways of doing this
        for matrix, power in zip(self.transforms, int_powers):
            for i in range(power):
                moment = moment @ matrix

        return moment

    def summation_form(self, moment) -> SummationSupercell:
        """ Convert to summation form """

        # TODO
