""" Base classes for PySpinW

This is an abstract outline, the actual implementations are in different files

"""

from abc import ABC, abstractmethod

import numpy as np
from ase.lattice import BravaisLattice


class MagneticStructure(ABC):
    """ Base class for representations of"""
    def __init__(self):
        pass


class Hamiltonian(ABC):
    """ Hamiltonian base class"""
    def __init__(self,
                 crystal_structure: BravaisLattice,
                 magnetic_structure: MagneticStructure):

        self.crystal_structure = crystal_structure
        self.magnetic_structure = magnetic_structure

    @abstractmethod
    def energies(self, q_vectors: np.ndarray):
        pass

class Sample(ABC):
    def __init__(self):
        pass