"""Selection of different Hamiltonians"""

from abc import ABC, abstractmethod

from pyspinw._base import MagneticStructure
from pyspinw.serialisation import SPWSerialisable
import numpy as np


# pylint: disable=R0903

class Hamiltonian(ABC, SPWSerialisable):
    """Hamiltonian base class"""

    def __init__(self,
                 magnetic_structure: MagneticStructure,
                 couplings: list[Coupling]):

        self.magnetic_structure = magnetic_structure
        self.couplings = couplings

    def get_structure(self):
        """Get the structural properties of the Hamiltonian needed for calculation."""
        raise NotImplementedError
        # sites = self.magnetic_structure.get_sites()
        # symmetry_stuff = self.magnetic_structure.get_symmetry()
        # etc....
        # return (sites, symmetry_stuff, etc...)

    @abstractmethod
    def energies(self, q_vectors: np.ndarray):
        """Calculate the energy levels of the system for the given q-vectors."""
        # this will have a specific implementation for each Hamiltonian
        raise NotImplementedError

    def serialise(self) -> dict:
        return {"magnetic_structure": self.magnetic_structure.serialise(),
                "couplings": [coupling.serialise() for coupling in self.couplings]}


class GeneralHamiltonian(Hamiltonian):
    """Hamiltonian of the form: H = Σ s.P.s + Q.s

    More explicitly, the equation is ... TODO

    """


class BiQuadratic(Hamiltonian):
    """Biquadratic Hamiltonian of the form H = - Σ A (s.s)^2

    More explicitly, the equation is ... TODO
    """
