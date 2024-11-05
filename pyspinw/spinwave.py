""" Direct python port of spinwave.m (version used saved in references) """

from dataclasses import dataclass
import numpy as np

import logging

logger = logging.getLogger("spinw")

@dataclass
class MagneticStructure:
    extended_magnetic_unit_cell_size: np.ndarray # 3 dimensional?
    propagation_vector_in_extended_cell_units: np.ndarray # 3 vector

    @property
    def propagation_vector(self):
        return self.propagation_vector_in_extended_cell_units * \
            self.extended_magnetic_unit_cell_size


    def is_commensurate(self, tol=1e-5) -> bool:
        # 394-405
        return np.any(np.abs(self.propagation_vector%1) > tol) # fractional part > tolerance

    def is_helical(self, tol=1e-5) -> bool:
        return np.any(np.abs(self.propagation_vector%1) > tol)

    def show_warnings(self):
        pass

@dataclass
class Lattice:
    """ Temporary object representing lattice details"""
    unit_cell_inverse_basis: np.ndarray

def spinwave(q_vectors: np.ndarray,
             magnetic_structure: MagneticStructure,
             lattice: Lattice):
    """
    spinwave.m port

    :param q_vectors: 3 x n_q vector of q points to evaluate at
    """

    # 279-284 - hkl is just an input as-is, call this q_vectors

    # 286-290 - warnings TODO: deal with warnings

    # 295-320 - ignore symbolic stuff

    # 322-327 - help, ignore this, python has its own things for this

    # 329-351 - title0 field: not sure we need this right now, or here
    #           inpForm: not sure we need this either

    # 353-365 - fastmode is just some other "usability" feature - don't need it

    # 367-384 - more guff

    # 387-405 - magnetic_structure input - in MagneticStructure
    #           get size of extended cell from structure, not sure we need to dereference
    #           warnings moved to magnetic structure class # TODO

    # 412-413 - Q vectors in lattice coordinates
    q_vectors_local = lattice.unit_cell_inverse_basis @ q_vectors

    # 415-416 - vectors in inverse angstroms - TODO: find out why and if needed

    # 418-420 - is helical, moved to MagneticStructure