""" Optimisation of propagation vectors """
from pyspinw.symmetry.supercell import RotationSupercell


class PropagationVectorOptimisation:
    def __init__(self, hamiltonian):

        if not isinstance(hamiltonian.structure.supercell, RotationSupercell):
            raise TypeError("We can only optimise the propagation vector in incommensurate structures "
                            "as they are required to take continuous values")

