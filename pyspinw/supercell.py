from fractions import Fraction
from math import lcm
from typing import Any

import numpy as np
from pydantic import BaseModel, root_validator, model_validator, field_validator

from pyspinw.tolerances import tolerances
from pyspinw.gui.cell_offsets import CellOffset

class PropagationVector(BaseModel):
    i: Fraction | float
    j: Fraction | float
    k: Fraction | float

    _vector: np.ndarray | None = None

    @property
    def vector(self):
        return self._vector

    def dot(self, cell_offset: CellOffset):
        return np.dot(self.vector, cell_offset.vector)

    def model_post_init(self, context: Any, /) -> None:
        self._vector = np.array([self.i, self.j, self.k], dtype=float)


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

    i: Fraction
    j: Fraction
    k: Fraction


    def model_post_init(self, context: Any, /) -> None:
        self._vector = np.array([self.i, self.j, self.k], dtype=float)

    def __init__(self, i: Fraction | float, j: Fraction | float, k: Fraction | float):
        super().__init__(i=_coerce_numeric_input(i),
                         j=_coerce_numeric_input(j),
                         k=_coerce_numeric_input(k))

class IncommensuratePropagationVector(PropagationVector):
    """ Propagation vector with non-integral values"""

    i: float
    j: float
    k: float

class CommensurateSupercell(BaseModel):

    propagation_vectors: list[CommensuratePropagationVector]


    model_config = {
        "arbitrary_types_allowed": True,
        "json_encoders": {
            np.ndarray: lambda v: v.tolist()
        }
    }

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
        raise NotImplementedError("summation_form not implemented in base class")

class SummationSupercell(CommensurateSupercell):
    """ Supercell with moments defined by

    m = sum(m_j exp(2 pi i r.k_j)).real
    """

    moments: list[np.ndarray] # TODO: not good for serialisation

    @field_validator("moments")
    def check_numpy_array(cls, list_of_arrays):
        for moment in list_of_arrays:
            if not isinstance(moment, np.ndarray):
                raise TypeError("moments must contain numpy ndarrays")
            if moment.shape != (3,):
                raise ValueError("moments must be of shape (3,)")

        return list_of_arrays

    @model_validator(mode="after")
    def check_lengths(self):
        if len(self.moments) != len(self.propagation_vectors):
            raise ValueError('the lists: moments and propagation_vector must have the same length')
        return self


    def evaluate(self, cell_offset: CellOffset, moment: np.ndarray):
        return np.sum(component.evaluate(cell_offset) for component in self.components).real

    def summation_form(self, moment) -> "SummationSupercell":
        return self

class TransformationSupercell(CommensurateSupercell):
    """ Supercell with moments defined by the following equation:

    m = prod(M_i^{r.k}) m0

    where m0 is the existing moment at a given site
    """

    transforms: list[np.ndarray]

    @field_validator("transforms")
    def check_numpy_array(cls, list_of_arrays):
        for moment in list_of_arrays:
            if not isinstance(moment, np.ndarray):
                raise TypeError("transforms must contain numpy ndarrays")
            if moment.shape != (3,3):
                raise ValueError("transforms must be of shape (3,3)")

        return list_of_arrays

    @model_validator(mode="after")
    def check_lengths(self):
        if len(self.transforms) != len(self.propagation_vectors):
            raise ValueError('the lists: transforms and propagation_vectors must have the same length')
        return self



    def evaluate(self, cell_offset: CellOffset, moment: np.ndarray):

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
