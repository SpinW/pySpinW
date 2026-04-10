"""pySpinW: Spin Hamiltonian calculations in Python

SpinW is a Python library that can plot and numerically simulate
magnetic structures and excitations
"""


from pyspinw.site import LatticeSite
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.structures import Structure
from pyspinw.anisotropy import Anisotropy, AxisMagnitudeAnisotropy
from pyspinw.exchange import HeisenbergExchange, DiagonalExchange, XYExchange, IsingExchange, DMExchange
from pyspinw.hamiltonian import Hamiltonian

from pyspinw.exchangegroup import InPlaneFilter, InDirectionFilter, BiDirectionFilter

from pyspinw.interface import (
    generate_sites, propagation_vectors, rotation_supercell, helical_supercell, summation_supercell, spacegroup,
    filter, generate_exchanges, axis_anisotropies, matrix_anisotropies, generate_helical_structure)

from pyspinw.symmetry.supercell import (
    SummationSupercell, RotationSupercell, TransformationSupercell, TiledSupercell,
    PropagationVector, CommensuratePropagationVector)

from pyspinw.sample import SingleCrystal, Multidomain, CrystalDomain, Twin, Powder, ScalingMethod
from pyspinw.path import Path, Path1D

from pyspinw.windows_parallelisation import set_up_windows_python_parallelisation

from pyspinw.calculations.spherical_integration import SphericalPointGeneratorType
from pyspinw.gui.viewer import show_hamiltonian as view

from pyspinw.units import CoordsUnits, IntensityUnits

from pyspinw.demo import demo_viewer, demo_chains
from pyspinw.demo import run_demos as demos

# TODO, add viewer and fitting things

__all__ = [
    # site
    "LatticeSite",

    # anisotropy
    "Anisotropy",
    "AxisMagnitudeAnisotropy",

    # exchanges
    "HeisenbergExchange",
    "DiagonalExchange",
    "XYExchange",
    "IsingExchange",
    "DMExchange",

    # filters
    "InPlaneFilter",
    "InDirectionFilter",
    "BiDirectionFilter",

    # structures / hamiltonian
    "Structure",
    "Hamiltonian",

    # symmetry
    "UnitCell",

    # units
    "CoordsUnits",
    "IntensityUnits",

    # interface functions
    "propagation_vectors",
    "rotation_supercell",
    "helical_supercell",
    "summation_supercell",
    "spacegroup",
    "filter",
    "axis_anisotropies",
    "matrix_anisotropies",
    "generate_sites",
    "generate_exchanges",
    "generate_helical_structure",

    # supercells
    "SummationSupercell",
    "RotationSupercell",
    "TransformationSupercell",
    "TiledSupercell",
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

    # viewer
    "view",

    # demo
    "demos",
    "demo_viewer",
    "demo_chains",
]
