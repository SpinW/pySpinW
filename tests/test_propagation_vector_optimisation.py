""" Tests for optimisation of propagation vectors"""
import pytest

from pyspinw.coupling import HeisenbergCoupling
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.site import LatticeSite
from pyspinw.structures import Structure
from pyspinw.symmetry.supercell import PropagationVector, RotationSupercell
from pyspinw.symmetry.unitcell import UnitCell

def test_end_to_end_matlab_test():
    """ Test from the MATLAB version

    function test_kbase(testCase, kbase_opts)
        % See https://doi.org/10.1103/PhysRevB.59.14367
        swobj = spinw();
        swobj.genlattice('lat_const', [3 3 8])
        swobj.addatom('r',[0; 0; 0],'S',1)
        swobj.gencoupling();
        J1 = 1.2;
        J2 = 1.0;
        swobj.addmatrix('label', 'J1', 'value', J1);
        swobj.addmatrix('label', 'J2', 'value', J2);
        swobj.addcoupling('mat', 'J1', 'bond', 2, 'subidx', 2);
        swobj.addcoupling('mat', 'J2', 'bond', 1);
        % Use rng seed for reproducible results
        swobj.optmagk('kbase', kbase_opts, 'seed', 1);

        expected_k = acos(-J2/(2*J1))/(2*pi);
        rel_tol = 1e-5;
        if abs(expected_k - swobj.mag_str.k(1)) > rel_tol*expected_k
            % If this k doesn't match, try 1-k
            expected_k = 1 - expected_k; % k and 1-k are degenerate
        end
        expected_mag_str = testCase.default_mag_str;
        expected_mag_str.k = [expected_k; expected_k; 0];
        testCase.verify_val(swobj.mag_str, expected_mag_str, ...
                            'rel_tol', 1e-5);
    end
    """

    j1 = 1.2
    j2 = 1.0


    unit_cell = UnitCell(3, 3, 8)
    s = LatticeSite(0,0,0,0,0,1)

    # Supercell to match the basis of [1 0; 0 1; 0 0]
    pv = PropagationVector(1, 1, 1/2)
    supercell = RotationSupercell([0,1,0], pv)

    structure = Structure([s], unit_cell, supercell=supercell)

    couplings = [HeisenbergCoupling(s, s, j=j1, name="J1", cell_offset=(1,0,0)),
                 HeisenbergCoupling(s, s, j=j2, name="J2", cell_offset=(0,1,0))]

    hamiltonian = Hamiltonian(structure, couplings)