"""pySpinW: Spin Hamiltonian calculations in Python

SpinW is a Python library that can plot and numerically simulate
magnetic structures and excitations
"""


from pyspinw.site import LatticeSite
from pyspinw.anisotropy import Anisotropy, AxisMagnitudeAnisotropy
from pyspinw.coupling import HeisenbergCoupling, DiagonalCoupling, XYCoupling, IsingCoupling, DMCoupling

from pyspinw.structures import Structure
from pyspinw.hamiltonian import Hamiltonian

from pyspinw.symmetry.unitcell import UnitCell

from pyspinw.interface import (
    sites, propagation_vectors, rotation_supercell, summation_supercell, spacegroup,
    filter, generate_exchanges, axis_anisotropies, matrix_anisotropies)

from pyspinw.symmetry.supercell import SummationSupercell, RotationSupercell, TransformationSupercell, TrivialSupercell

from pyspinw.sample import SingleCrystal, Multidomain, CrystalDomain, Twin, Powder, ScalingMethod
from pyspinw.path import Path, Path1D

from windows_parallelisation import set_up_windows_python_parallelisation

# TODO, add viewer and fitting things

__all__ = [
    # site
    "LatticeSite",

    # anisotropy
    "Anisotropy",
    "AxisMagnitudeAnisotropy",

    # coupling
    "HeisenbergCoupling",
    "DiagonalCoupling",
    "XYCoupling",
    "IsingCoupling",
    "DMCoupling",

    # structures / hamiltonian
    "Structure",
    "Hamiltonian",

    # symmetry
    "UnitCell",

    # interface functions
    "sites",
    "propagation_vectors",
    "rotation_supercell",
    "summation_supercell",
    "spacegroup",
    "filter",
    "generate_exchanges",
    "axis_anisotropies",
    "matrix_anisotropies",

    # supercells
    "SummationSupercell",
    "RotationSupercell",
    "TransformationSupercell",
    "TrivialSupercell",

    # sample
    "SingleCrystal",
    "Multidomain",
    "CrystalDomain",
    "Twin",
    "Powder",
    "ScalingMethod",

    # path
    "Path",
    "Path1D",

    # windows parallelisation
    "set_up_windows_python_parallelisation",
]