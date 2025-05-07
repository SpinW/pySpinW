"""Base classes for PySpinW

This is an abstract outline, the actual implementations are in different files

"""

# pylint: disable=R0903

from abc import ABC, abstractmethod

import numpy as np
from ase.lattice import BravaisLattice


class MagneticStructure(ABC):
    """Base class for representations of the Magnetic Structures"""

    def __init__(self):
        pass


class Hamiltonian(ABC):
    """Hamiltonian base class"""

    def __init__(self,
                 crystal_structure: BravaisLattice,
                 magnetic_structure: MagneticStructure):

        self.crystal_structure = crystal_structure
        self.magnetic_structure = magnetic_structure

    @abstractmethod
    def energies(self, q_vectors: np.ndarray):
        """Get the energy levels of the system at the given q vectors"""

class Sample(ABC):
    """Representation of the macrostructure of a sample used in an experiment (Twin, Powder etc)"""

    def __init__(self, hamiltonian: Hamiltonian):
        self.hamiltonian = hamiltonian

Identifier = str # temporary choice for now


class Coupling:
    """Coupling between different sites"""

    def __init__(self, site_1: Identifier, site_2: Identifier):
        self._site_1 = site_1
        self._site_2 = site_2
        self._coupling_matrix = None

    @property
    def coupling_matrix(self) -> np.ndarray:
        """The coupling matrix for this coupling

        i.e. if H is the energy contribution for this coupling, S is the spin state, and
        M is the coupling matrix, we have

        H = S^T M S
        """
        if self._coupling_matrix is None:
            raise ValueError("Coupling matrix not initialised - this shouldn't happen")
        else:
            return self._coupling_matrix


class Anisotropy:
    """Defines the anisotropy at a given site"""

    def __init__(self, site: Identifier):
        self._site = site
        self._anisotropy_matrix = None

    @property
    def anisotropy_matrix(self) -> np.ndarray:
        """Matrix spefifying the anisotropy - `A` term in the Hamiltonian"""
        if self._anisotropy_matrix is None:
            raise ValueError("Anisotropy matrix not initialised - this shouldn't happen")
        else:
            return self._anisotropy_matrix


class Data:
    """Placeholder"""

class Experiment:
    """The setup of a neutron experiment (base class)"""

    def __init__(self, sample: Sample, data: Data | None = None):
        pass
