"""Different Kinds of samples"""

from abc import ABC, abstractmethod

from pyspinw.hamiltonian import Hamiltonian

# pylint: disable=R0903

class Sample(ABC):
    """Representation of the macrostructure of a sample used in an experiment (Twin, Powder etc)"""

    def __init__(self, hamiltonian: Hamiltonian):
        self.hamiltonian = hamiltonian

    @abstractmethod
    def generate_q(self, input_q: list[np.ndarray], n_q: int, resolution, field):
        """Generate a sample of q-vectors for a given magnetic field and resolution."""
        raise NotImplementedError

class SingleCrystal(Sample):
    """Specifies a single crystal sample"""


class Multidomain(Sample):
    """Sample consisting of multiple domains"""


class Twin(Multidomain):
    """Specify a twinned crystal.

    Special case of multidomain
    """


class Powder(Sample):
    """Sample is a powder"""
