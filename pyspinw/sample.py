"""Different Kinds of samples"""

from abc import ABC, abstractmethod

import numpy as np
from numpy._typing import ArrayLike

from docs.developers.workflow_example_prototype import hamiltonian
from pyspinw.checks import check_sizes
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.path import Path, Path1D
from pyspinw.tolerances import tolerances


# pylint: disable=R0903

class Sample(ABC):
    """Representation of the macrostructure of a sample used in an experiment (Twin, Powder etc)"""

    def __init__(self, hamiltonian: Hamiltonian):
        self.hamiltonian = hamiltonian

    @abstractmethod
    def generate_q(self, input_q: list[np.ndarray], n_q: int, resolution, field):
        """Generate a sample of q-vectors for a given magnetic field and resolution."""
        raise NotImplementedError()

    @abstractmethod
    def _energies(self, q_values: np.ndarray, use_rust: bool=True):
        """ Get the energies as a function of the q samples"""
        return self.hamiltonian.energies(q_values, use_rust=use_rust)

    @abstractmethod
    def _q_energy_intensity_samples(self) -> np.ndarray:
        """Returns a list of sampled intensities along with the q value and energy it corresponds to"""
        raise NotImplementedError()

class Sample3D(Sample):
    """ Sample where the direction of q matters"""
    @abstractmethod
    def energies(self,
                 path: Path,
                 field: ArrayLike | None = None,
                 use_rust: bool=True) -> list[np.ndarray]:

        """ Get the energies along a specified path """


class Sample1D(Sample):
    """ Sample where the direction of q does not matter"""
    def energies(self, path: Path1D, field: ArrayLike | None = None, use_rust: bool=True):
        pass


class SingleCrystal(Sample3D):
    """Specifies a single crystal sample"""
    def __init__(self, hamiltonian: Hamiltonian):

        super().__init__(hamiltonian=hamiltonian)

    def energies(self,
                 path: Path,
                 field: ArrayLike | None = None,
                 use_rust: bool=True):

        return self._energies(path.q_points())


class CrystalDomain:
    """ Orientation and amount of crystalline domain (e.g. one of a twin, or a domain in the traditional sense) """
    @check_sizes(transformation=(3,3), force_numpy=True)
    def __init__(self, transformation: ArrayLike, weighting: float):

        # check the transformation is valid
        if not np.allclose(transformation @ transformation.T, np.eye(3), atol=tolerances.IS_ZERO_TOL):
            raise ValueError(f"Expected crystal domain transformation to be orthonormal, got {transformation}")

        if weighting < 0:
            raise ValueError(f"Weightings must be positive, got {weighting}")

        self.transformation = transformation
        self.weighting = weighting

CrystalDomainLike = CrystalDomain | tuple[ArrayLike, float]


def _domain_like_to_domain(domain_like: CrystalDomainLike) -> CrystalDomain:
    """ Converts CrystalDomainLike to CrystalDomain """

    if isinstance(domain_like, CrystalDomain):
        return domain_like
    else:
        return CrystalDomain(*domain_like)


class Multidomain(Sample3D):
    """Sample consisting of multiple domains"""

    def __init__(self,
                 hamiltonian: Hamiltonian,
                 domains: list[CrystalDomainLike]):

        self._domains = [_domain_like_to_domain(domain) for domain in domains]

        self._transformations = [domain.transformation for domain in self._domains]

        weights = [domain.weighting for domain in self._domains]
        total_weight = sum(weights)
        self.weights = [weight/total_weight for weight in weights]

        super().__init__(hamiltonian=hamiltonian)

    @check_sizes(field=(3,), allow_nones=True)
    def energies(self,
                 path: Path,
                 field: ArrayLike | None = None,
                 use_rust: bool=True) -> list[np.ndarray]:

        q_vectors = path.q_points()

        output = []
        for transformation in self._transformations:

            transformed_q = transformation @ q_vectors
            transformed_field = None if field is None else transformation @ field

            output += self.hamiltonian.energies(transformed_q, field=transformed_field, use_rust=use_rust)

        return output



class Twin(Multidomain):
    """Specify a twinned crystal.

    Special case of multidomain
    """

    def __init__(self,
                 hamiltonian: Hamiltonian,
                 relative_rotation: ArrayLike,
                 second_twin_fraction: float):

        if not (0 <= second_twin_fraction <= 1):
            raise ValueError("Expected weighting to be between 0 and 1")

        first_twin_fraction = 1 - second_twin_fraction

        super().__init__(
            hamiltonian=hamiltonian,
            domains=[CrystalDomain(np.eye(3), first_twin_fraction),
                     CrystalDomain(relative_rotation, second_twin_fraction)])

class Powder(Sample1D):
    """Sample is a powder"""


class TwoMagnon(Sample3D):
    """ Two magnon excitations in a single crystal"""
