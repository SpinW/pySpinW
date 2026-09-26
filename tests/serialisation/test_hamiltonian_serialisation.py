""" Tests for serialisation of the Hamiltonian class """
import numpy as np
import pytest

from pyspinw import *

def test_hamiltonian_serialisation():

    sites = [
        LatticeSite(0.1, 0.2, 0.3, 0, 0, 1),
        LatticeSite(0.5, 0.5, 0, 1, 0, 0)
    ]

    sg = spacegroup("P4mmm")
    unit_cell = UnitCell(1,1,1)
    supercell = TiledSupercell(scaling=(2,3,5))

    structure = Structure(sites, unit_cell=unit_cell, supercell=supercell, spacegroup=sg)

    anisotropies = [AxisMagnitudeAnisotropy(sites[0], direction=(1,0,0), a=-1),
                    Anisotropy(sites[1], anisotropy_matrix=np.eye(3))]

    exchanges = []

