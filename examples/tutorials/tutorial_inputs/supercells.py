""" Tutorial about creating different kinds of supercells """

## title: Supercells

# Specifying magnetic structure at scales beyond the crystallographic unit cell is done using the `Supercell`
# classes. The supercell allows the user to specify a pattern of spins across different unit cells, but leaves
# everything else (atom positions, exchanges, anisotropies) unchanged.

## subtitle: Commensurate and Incommensurate Supercells

# There are two categories of supercell available, which we call "commensurate" and "incommensurate", this
# distinction is made because of the different strategies that pySpinW uses for calculating different kinds
# of system. The two strategies are this
#
#   (1) internally, spinW can explicitly expand the supercell, creating copies of the original cell,
#       and increasing the number of spins in the calculation by the number of cells in the supercell.
#       This can only be done when the supercell is a rational multiple of the unit cell, or in other words,
#       when the propagation vectors are commensurate. In comparison with the alternative strategy,
#       this can be relatively slow as we increase the size of the system, and it can introduce
#       "ghost modes", where there are magnon modes with zero intensity. However, this kind of
#       supercell is more flexible in other ways.
#   (2) spinW can use the method from [Toth and Lake](https://iopscience.iop.org/article/10.1088/0953-8984/27/16/166002)
#       to *implicitly* calculate magnon energies and scattering intensities for *a single*, potentially
#       incommensurate, propagation vector.
#       This can be cheaper and doesn't have (as many) ghost modes.
#       It can be used for improving calculations with a single commensurate propagation vector, as long as the
#       magnetic structure can be described in this way.
#
# There is, therefore, a choice that one must make about which is most appropriate for any given system.

## subtitle: Single and Multiple Component Spin Definitions

# The is another kind of distinction one can make between types of supercell.
#
#
#
#
# The remainder of this tutorial runs through the different supercell options that are available in spinW.
#
# As usual, we'll do a `*` import of the main pySpinW features

from pyspinw import *

## subtitle: Tiled Supercell (Commensurate)

# This is the simplest supercell and just makes exact copies of the original cell. It doesn't have much
# real-world value, but can be useful for debugging, and for showing how supercells work.
#
# As an example, we will make a structure

structure = Structure(
    [LatticeSite(0.5, 0.5, 0.5, 0, 0, 1)],
    UnitCell(1, 1, 1),
    spacegroup("p1"),
    supercell=TiledSupercell(scaling=(3,6,1)))

view(structure)
