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
                k: Fraction | float,
                phase: float = 0.0):

        self.i = i
        self.j = j
        self.k = k
        self.phase = phase

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
            "phase": self.phase,
            "assured_commensurate": False
        }

    @staticmethod
    @expects_keys("i,j,k,phase")
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        i = deserialise_fraction_or_builtin(json["i"])
        j = deserialise_fraction_or_builtin(json["j"])
        k = deserialise_fraction_or_builtin(json["k"])
        phase = json["phase"]

        if json["assured_commensurate"]:
            return CommensuratePropagationVector(i, j, k, phase=phase)
        else:
            return PropagationVector(i, j, k, phase=phase)

    def __eq__(self, other):
        return self.i == other.i and self.j == other.j and self.k == other.k

    def uncorrected_phase_position(self, site: LatticeSite):
        """ Get the position of a site along the propagation vector, from (0,0,0), ignoring the phase correction"""

        return 2*np.pi*(self.vector * site.ijk)

    def corrected_phase_position(self, site: LatticeSite):
        """ Get the position of a site along the propagation vector, from (0,0,0), including the phase correction"""

        return 2 * np.pi * (self.vector * site.ijk) + self.phase

class CommensuratePropagationVector(PropagationVector):
    """ Propagation vector with integer values"""

    def __init__(self,
                 i: Fraction | float,
                 j: Fraction | float,
                 k: Fraction | float,
                 phase: float = 0.0):

        i = _coerce_numeric_input(i)
        j = _coerce_numeric_input(j)
        k = _coerce_numeric_input(k)

        super().__init__(i, j, k, phase=phase)

        self._vector = np.array([self.i, self.j, self.k], dtype=float)


    def _serialise(self, context: SPWSerialisationContext) -> dict:

        # Just use the parent class serialiser, but set the 'assured_commensurate' flag
        serialised = super()._serialise(context)
        serialised["assured_commensurate"] = True
        return serialised

    def __eq__(self, other):
        if isinstance(other, CommensuratePropagationVector):
            return self.i == other.i and self.j == other.j and self.k == other.k

        return False

class SupercellTransformation(ABC, SPWSerialisable):
    """ Base class for supercell transformations """

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
    """ Transformation that does nothing """

    transformation_name = "identity"

    def _serialise_transform(self) -> dict:
       return {}

    @staticmethod
    def _deserialise_transform(json) -> "SupercellTransformation":
        return IdentityTransform()

    def apply(self, moment: np.ndarray, propagation_vector: CommensuratePropagationVector, cell_offset: CellOffset):
        """ Apply this transformation to a given moment (as it is the identity it does nothing to it)"""
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
        """ Apply this transformation (rotation) to a moment in a specified cell """
        angle = 2 * np.pi * propagation_vector.dot(cell_offset) + propagation_vector.phase
        return moment @ rotation_matrix(angle, self._axis)

transform_types = {cls.transformation_name: cls for cls in [IdentityTransform, RotationTransform]}


class Supercell(ABC, SPWSerialisable):
    """ Base class for different supercells """

    serialisation_name = "supercell"
    supercell_name = "<base supercell>"

    def __init__(self, scaling: tuple[int, int, int] | None = None):

        if scaling is None:
            self._scaling = (1,1,1)
        else:
            self._scaling = scaling

        if any([scaling_component < 1 for scaling_component in self._scaling]):
            raise ValueError("Scaling components should be >= 1")

    @abstractmethod
    def moment(self, site: LatticeSite, cell_offset: CellOffset):
        """ Evaluate the moment of a lattice site in this unit cell """

    def cell_position_and_moment(self, site: LatticeSite, cell_offset: CellOffset):
        """ Position within the unit cell, and moment"""
        return np.concatenate((site.ijk, self.moment(site, cell_offset)))

    def supercell_position_and_moment(self, site: LatticeSite, cell_offset: CellOffset):
        """ Position within the supercell and moment """
        return np.concatenate((site.ijk + cell_offset.vector, self.moment(site, cell_offset)))

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
            supercell_type = json["type"]
            further_data = json["data"]
            return supercell_types[supercell_type]._deserialise_supercell(further_data, scale, context)

        except KeyError as ke:
            expected_names = ", ".join([f"'{key}'" for key in supercell_types])
            raise SPWSerialisationError(f"Expected transform type to be one of {expected_names}") from ke


class TrivialSupercell(Supercell):
    """ Trivial supercell, just a single unit cell """

    supercell_name = "trivial"

    def moment(self, site: LatticeSite, cell_offset: CellOffset):
        """ Get the moment for a site with an offset within the supercell """
        return site.base_moment

    def cell_size(self) -> tuple[int, int, int]:
        """ How big is this supercell """
        return self._scaling

    def summation_form(self) -> "Supercell":
        """ Get this supercell in summation form """
        return self

    def _serialise_supercell(self, context: SPWSerialisationContext):
        return {}

    @staticmethod
    def _deserialise_supercell(json, scale, context: SPWDeserialisationContext):
        return TrivialSupercell(scale)


class CommensurateSupercell(Supercell):
    """ Base class for a commensurate supercell"""

    supercell_name = "<commensurate supercell base>"

    def __init__(self,
                 propagation_vectors: list[CommensuratePropagationVector],
                 scaling: tuple[int, int, int] | None = None):

        self._propagation_vectors = propagation_vectors
        super().__init__(scaling)


    def _transform_evaluate(self, cell_offset: CellOffset, moment: np.ndarray):
        """ Apply the transformation to a given moment """
        raise NotImplementedError("evaluate not implemented in base class")

    def moment(self, site: LatticeSite, cell_offset: CellOffset):
        """ Get the moment for a site in a specified cell according to the supercell"""
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

    supercell_name = "transformation"

    def __init__(self,
                 transforms: list[tuple[CommensuratePropagationVector, SupercellTransformation | None]],
                 scaling=(1, 1, 1)):

        # TODO: Provide a nicer interface for this maybe
        self._transforms = [(vector, IdentityTransform() if transform is None else transform)
                            for vector, transform in transforms]

        propagation_vectors = [vector for vector, _ in transforms]

        super().__init__(propagation_vectors, scaling)

    def _transform_evaluate(self, cell_offset: CellOffset, moment: np.ndarray):
        for vector, transform in self._transforms:
            moment = transform.apply(moment=moment, propagation_vector=vector, cell_offset=cell_offset)

        return moment

    def summation_form(self) -> "Supercell":
        """ Convert into summation form """
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

        return TransformationSupercell(transforms, scale)


class SummationSupercell(CommensurateSupercell):
    """ Supercell with moments defined by

    m = sum(m_j exp(2 pi i r.k_j)).real
    """

    supercell_name = "summation"

    def __init__(self,
                 propagation_vectors: list[CommensuratePropagationVector],
                 scaling: tuple[int, int, int] | None = None):

        super().__init__(propagation_vectors, scaling)

    def moment(self, site: LatticeSite, cell_offset: CellOffset):
        """ Calculate moment at a given cell offset"""
        moment = np.zeros(3, dtype=complex)
        for component, propagation_vector in zip(site.moment_data, self._propagation_vectors):
            # TODO: check
            moment += component * np.exp(2j * np.pi * propagation_vector.dot(cell_offset))

        return moment.real

    def summation_form(self) -> "SummationSupercell":
        """ Convert to summation form (it is already in this form, but not all Supercells are)"""
        return self

    def _serialise_supercell(self, context: SPWSerialisationContext):
        return {"vectors": [vector._serialise(context) for vector in self._propagation_vectors]}

    @staticmethod
    @expects_keys("vectors")
    def _deserialise_supercell(json, scale: tuple[int, int, int], context: SPWDeserialisationContext):
        vectors = [PropagationVector._deserialise(data, context) for data in json["vectors"]]
        return SummationSupercell(vectors, scale)

supercell_types = {cls.supercell_name: cls
                   for cls in [TrivialSupercell, TransformationSupercell, SummationSupercell]}
