"""pySpinW: Spin Hamiltonian calculations in Python

SpinW is a Python library that can plot and numerically simulate
magnetic structures and excitations
"""


from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structures import Structure
from pyspinw.anisotropy import Anisotropy, AxisMagnitudeAnisotropy
from pyspinw.coupling import HeisenbergCoupling, DiagonalCoupling, XYCoupling, IsingCoupling, DMCoupling
from pyspinw.hamiltonian import Hamiltonian

from pyspinw.couplinggroup import InPlaneFilter, InDirectionFilter, BiDirectionFilter

from pyspinw.interface import (
    generate_sites, propagation_vectors, rotation_supercell, summation_supercell, spacegroup,
    filter, generate_exchanges, axis_anisotropies, matrix_anisotropies)

from pyspinw.symmetry.supercell import (
    SummationSupercell, RotationSupercell, TransformationSupercell, TrivialSupercell,
    PropagationVector, CommensuratePropagationVector)

from pyspinw.sample import SingleCrystal, Multidomain, CrystalDomain, Twin, Powder, ScalingMethod
from pyspinw.path import Path, Path1D

from pyspinw.windows_parallelisation import set_up_windows_python_parallelisation

from pyspinw.calculations.spherical_integration import SphericalPointGeneratorType

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

    # coupling groups
    "InPlaneFilter",
    "InDirectionFilter",
    "BiDirectionFilter",

    # structures / hamiltonian
    "Structure",
    "Hamiltonian",

    # symmetry
    "UnitCell",

    # interface functions
    "generate_sites",
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
    "PropagationVector",
    "CommensuratePropagationVector",

    # sample
    "SingleCrystal",
    "Multidomain",
    "CrystalDomain",
    "Twin",
    "Powder",
    "ScalingMethod",

    # Point generators
    "SphericalPointGeneratorType",

    # path
    "Path",
    "Path1D",

    # windows parallelisation
    "set_up_windows_python_parallelisation",
]
