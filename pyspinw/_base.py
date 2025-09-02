"""Base classes for PySpinW

This is an abstract outline, the actual implementations are in different files

"""

# pylint: disable=R0903

from abc import ABC, abstractmethod
from typing import ClassVar

import numpy as np
from ase.lattice import BravaisLattice
from pydantic import BaseModel

from pyspinw.checks import check_sizes
from pyspinw.gui.cell_offsets import CellOffset
from pyspinw.serialisation import SPWSerialisable, SPWSerialisationContext, numpy_serialise, numpy_deserialise
from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell


class MagneticStructure(SPWSerialisable):
    """Base class for representations of the Magnetic Structures"""

    def __init__(self):
        pass

    def serialise(self):
        return {}

    def deserialise(data: dict):
        return MagneticStructure()

class Coupling(SPWSerialisable):
    """Coupling between different sites"""

    def __init__(self,
                 name: str,
                 site_1: LatticeSite,
                 site_2: LatticeSite,
                 cell_offset: CellOffset | None):

        self.name = name
        self.site_1 = site_1
        self.site_2 = site_2

        self.cell_offset = CellOffset(i=0,j=0,k=0) if cell_offset is None else cell_offset


    coupling_type: ClassVar[str] = "Base Coupling"
    parameters: ClassVar[list[str]] = []
    parameter_defaults: ClassVar[list[int]] = []
    short_string: ClassVar[str] = "X"

    @property
    def coupling_matrix(self) -> np.ndarray:
        """The coupling matrix for this coupling

        i.e. if H is the energy contribution for this coupling, S is the spin state, and
        M is the coupling matrix, we have

        H = S^T M S
        """
        raise NotImplementedError()

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



class Data:
    """Placeholder"""

    def __init__(self, data):
        self.data = data

    @property
    def q(self) -> np.ndarray:
        raise NotImplementedError
