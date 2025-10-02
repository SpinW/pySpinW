""" Supercell representations"""
from abc import ABC, abstractmethod
from fractions import Fraction
from math import lcm

import numpy as np

from pyspinw.checks import check_sizes
from pyspinw.serialisation import SPWSerialisable, SPWSerialisationContext, \
    serialise_fraction_or_builtin, deserialise_fraction_or_builtin, SPWDeserialisationContext, expects_keys, \
    SPWSerialisationError, numpy_serialise, numpy_deserialise
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

        # Just use the parent class serialiser, but set the 'assured_commensurate' flag
        serialised = super()._serialise(context)
        serialised["assured_commensurate"] = True
        return serialised

class SupercellTransformation(ABC, SPWSerialisable):

    serialisation_name = "supercell_transformation"
    transformation_name = "<transform base class>"

    @abstractmethod
    def apply(self, moment: np.ndarray, propagation_vector: CommensuratePropagationVector, cell_offset: CellOffset):
        """ Apply this transformation to a moment with a given offset """

    @abstractmethod
    def _serialise_transform(self) -> dict:
        """ Serialise this transform"""

    @staticmethod
    @abstractmethod
    def _deserialise_transform(json) -> "SupercellTransformation":
        """ Deserialise this transform"""


    def _serialise(self, context: SPWSerialisationContext) -> dict:
        return {
            "type": self.transformation_name,
            "data": self._serialise_transform()
        }

    @staticmethod
    @expects_keys("type, data")
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        try:
            return transform_types[json["type"]]._deserialise_transform(json["data"])

        except KeyError as ke:
            expected_names = ", ".join([f"'{key}'" for key in transform_types])
            raise SPWSerialisationError(f"Expected transform type to be one of {expected_names}") from ke


class IdentityTransform(SupercellTransformation):

    transformation_name = "identity"

    def _serialise_transform(self) -> dict:
       return {}

    @staticmethod
    def _deserialise_transform(json) -> "SupercellTransformation":
        return IdentityTransform()

    def apply(self, moment: np.ndarray, propagation_vector: CommensuratePropagationVector, cell_offset: CellOffset):
        return moment

class RotationTransform(SupercellTransformation):
    """ Rotation transformation in a supercell"""

    transformation_name = "rotation"

    @check_sizes(axis=(3,), force_numpy=True)
    def __init__(self, axis: np.ndarray | list[float]):
        self._axis = axis

    def _serialise_transform(self):
        return {"axis": numpy_serialise(self._axis)}

    @staticmethod
    @expects_keys("axis")
    def _deserialise_transform(json) -> "SupercellTransformation":
        return RotationTransform(axis=numpy_deserialise(json["axis"]))

    def apply(self, moment: np.ndarray, propagation_vector: CommensuratePropagationVector, cell_offset: CellOffset):
        angle = 2 * np.pi * propagation_vector.dot(cell_offset)
        return moment @ rotation_matrix(angle, self._axis)

transform_types = {cls.transformation_name: cls for cls in [IdentityTransform, RotationTransform]}


class Supercell(ABC, SPWSerialisable):

    serialisation_name = "supercell"
    supercell_name = "<base supercell>"

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

    def cells(self):
        """ Iterator for cells in supercell """
        a, b, c = self.cell_size()
        for i in range(a):
            for j in range(b):
                for k in range(c):
                    yield CellOffset(i,j,k)

    @abstractmethod
    def summation_form(self) -> "Supercell":
        """Get a summation type supercell"""

    @abstractmethod
    def _serialise_supercell(self, context: SPWSerialisationContext):
        """ Serialise this supercell """

    @staticmethod
    @abstractmethod
    def _deserialise_supercell(json, scale: tuple[int, int, int], context: SPWDeserialisationContext):
        """ Deserialise this supercell type """

    def _serialise(self, context: SPWSerialisationContext) -> dict:
        return {
            "type": self.supercell_name,
            "data": self._serialise_supercell(context),
            "s_i": self._scaling[0],
            "s_j": self._scaling[1],
            "s_k": self._scaling[2]
        }

    @staticmethod
    @expects_keys("type, data, s_i, s_j, s_k")
    def _deserialise(json: dict, context: SPWDeserialisationContext):

        try:
            scale = (json["s_i"], json["s_j"], json["s_k"])
            return supercell_types[json["type"]]._deserialise_supercell(json["data"], scale, context)

        except KeyError as ke:
            expected_names = ", ".join([f"'{key}'" for key in supercell_types])
            raise SPWSerialisationError(f"Expected transform type to be one of {expected_names}") from ke


class TrivialSupercell(Supercell):
    """ Trivial supercell, just a single unit cell """

    supercell_name = "trivial"

    def moment(self, site: LatticeSite, cell_offset: CellOffset):
        return site.base_moment

    def cell_size(self) -> tuple[int, int, int]:
        return self._scaling

    def summation_form(self) -> "Supercell":
        return self

    def _serialise_supercell(self, context: SPWSerialisationContext):
        return {}

    @staticmethod
    def _deserialise_supercell(json, scale, context: SPWDeserialisationContext):
        return TrivialSupercell(scale)


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
        i = lcm(*[vector.i.denominator for vector in self._propagation_vectors])
        j = lcm(*[vector.j.denominator for vector in self._propagation_vectors])
        k = lcm(*[vector.k.denominator for vector in self._propagation_vectors])

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

    def __init__(self, transforms: list[tuple[CommensuratePropagationVector, SupercellTransformation | None]], scale=(1,1,1)):
        # TODO: Provide a nicer interface for this maybe
        self._transforms = [(vector, IdentityTransform() if transform is None else transform)
                            for vector, transform in transforms]

        propagation_vectors = [vector for vector, _ in transforms]

        super().__init__(propagation_vectors, scale)

    def _transform_evaluate(self, cell_offset: CellOffset, moment: np.ndarray):
        for vector, transform in self._transforms:
            moment = transform.apply(moment=moment, propagation_vector=vector, cell_offset=cell_offset)

        return moment

    def summation_form(self) -> "Supercell":
        raise NotImplementedError("Not implemented yet")

    def _serialise_supercell(self, context: SPWSerialisationContext):
        return [{"vector": vector._serialise(context),
                 "transform": transform._serialise(context)}
                    for vector, transform in self._transforms]

    @staticmethod
    def _deserialise_supercell(json, scale, context):
        transforms = [(PropagationVector._deserialise(part["vector"], context),
                       SupercellTransformation._deserialise(part["transform"], context))
                      for part in json]

        return TransformationSupercell(scale, transforms)


class SummationSupercell(CommensurateSupercell):
    """ Supercell with moments defined by

    m = sum(m_j exp(2 pi i r.k_j)).real
    """

    def __init__(self, propagation_vectors):
        super().__init__(propagation_vectors)


    def moment(self, cell_offset: CellOffset, moment: np.ndarray):
        """ Calculate moment at a given cell offset"""
        return np.sum(component.evaluate(cell_offset) for component in self.components).real

    def summation_form(self) -> "SummationSupercell":
        """ Convert to summation form (it is already in this form, but not all Supercells are)"""
        return self

    def _serialise_supercell(self, context: SPWSerialisationContext):
        raise NotImplementedError("Serialisation not implemented") # TODO

    @staticmethod
    def _deserialise_supercell(json, scale: tuple[int, int, int], context: SPWDeserialisationContext):
        raise NotImplementedError("Serialisation not implemented") # TODO

supercell_types = {cls.supercell_name: cls
                   for cls in [TrivialSupercell, TransformationSupercell, SummationSupercell]}