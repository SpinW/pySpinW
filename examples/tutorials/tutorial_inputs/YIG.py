
## title: Yttrium Iron Garnet Example

# This tutorial uses a CIF file for specifying the crystal structure, [download it here](YIG_10K.cif)
## file: YIG_10K.cif

# Import pyspinw and load in the structure

from pyspinw import *

structure = load("YIG_10K.cif")

view(structure)
