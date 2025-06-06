"""Base classes for PySpinW

This is an abstract outline, the actual implementations are in different files

"""

# pylint: disable=R0903

from abc import ABC, abstractmethod
from typing import ClassVar

import numpy as np
from ase.lattice import BravaisLattice
from pydantic import BaseModel

from pyspinw.gui.cell_offsets import CellOffset
from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell


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


class Coupling(BaseModel):
    """Coupling between different sites"""

    name: str

    site_1: LatticeSite
    site_2: LatticeSite

    cell_offset: CellOffset

    coupling_type: ClassVar[str] = "Base Coupling"
    parameters: ClassVar[list[str]] = []
    parameter_defaults: ClassVar[list[int]] = []
    short_string: ClassVar[str] = "X"

    _coupling_matrix: np.ndarray | None = None


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

    @property
    def parameter_string(self) -> str:
        """ String representation of parameters """
        return ", ".join([f"{parameter}={self.__dict__[parameter]:.5g}" for parameter in self.parameters])

    @property
    def lattice_vector(self):
        """ Vector from site 1 to site 2 in lattice coordinates"""
        return self.cell_offset.as_tuple + self.site_2.ijk - self.site_1.ijk

    def vector(self, unit_cell: UnitCell):
        """ Vector from site 1 to site 2 in cartesian coordinates (requires a unit cell definition)"""
        return unit_cell.fractional_to_cartesian(self.lattice_vector)

    def distance(self, unit_cell: UnitCell):
        """ Distance between sites """
        return np.sqrt(np.sum(self.vector(unit_cell)))

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
