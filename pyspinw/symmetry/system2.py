""" Lattice systems

Classes and data for working with lattice systems
"""

class LatticeSystem:
    pass

class Trigonal(LatticeSystem):
    pass

class Monoclinic(LatticeSystem):
    def __init__(self, unique_axis = "b"):
        self.unique_axis = unique_axis

class Orthorhombic(LatticeSystem):
    pass

class Tetragonal(LatticeSystem):
    pass

class Cubic(LatticeSystem):
    pass