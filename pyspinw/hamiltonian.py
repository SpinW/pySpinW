"""Selection of different Hamiltonians"""

from abc import ABC, abstractmethod

from pyspinw.coupling import Coupling
from pyspinw.serialisation import SPWSerialisable
import numpy as np

from pyspinw.structures import Structure


# pylint: disable=R0903

class Hamiltonian(ABC, SPWSerialisable):
    """Hamiltonian base class"""

    serialisation_name = "hamiltonian"

    def __init__(self,
                 structure: Structure,
                 couplings: list[Coupling]):

        self.structure = structure
        self.couplings = couplings

    @abstractmethod
    def energies(self, q_vectors: np.ndarray):
        """Calculate the energy levels of the system for the given q-vectors."""



    def _serialise(self) -> dict:
        return {"magnetic_structure": self.structure._serialise(),
                "couplings": [coupling._serialise() for coupling in self.couplings]}

