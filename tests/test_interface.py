""" Tests for interface functions (doubles as system tests) """

# Tests compares with Matlab output, generated from: 
# https://gist.github.com/mducle/1026b0940f3f9cf4f722261e728885d7

import pytest
import numpy as np
import scipy.io
import os

from pyspinw.hamiltonian import Hamiltonian, egrid, omegasum
from pyspinw.interface import generate_exchanges, generate_structure, generate_helical_structure, filter, axis_anisotropies
from pyspinw.symmetry.unitcell import UnitCell

QVECS = np.array([[0.1, 0, 0], [0.2, 0, 0], [0.1, 0.1, 0], [0.1, 0.1, 0.1], [0, 0.1, 0.1], [0, 0, 0.1]])
REFDAT = scipy.io.loadmat(os.path.join(os.path.dirname(__file__), 'matlab_test_data.mat'))
INTHIGH = 100

# Monkey patch the parallelisation routines so calculations show up in coverage
from concurrent.futures import ThreadPoolExecutor
import sys
def get_Executor():
    return ThreadPoolExecutor
sys.modules['pyspinw.calculations.spinwave'].get_Executor = get_Executor

def _test_ref_data(test_name, energy, intensity, ignoreQ=None):
    eref, intref = omegasum(REFDAT[test_name][0][0].T, REFDAT[test_name][0][1].T, is_series=False)
    # Remove very large intensities as this diverges near zone centre
    intensity[np.where(intensity > INTHIGH)] = np.nan 
    intref[np.where(intref > INTHIGH)] = np.nan 
    if ignoreQ is not None:
        idx = [i for i in range(energy.shape[0]) if i not in ignoreQ]
        energy, intensity, eref, intref = (energy[idx,:], intensity[idx,:], eref[idx,:], intref[idx,:])
    np.testing.assert_allclose(energy, eref, atol=1e-5, rtol=0)
    np.testing.assert_allclose(intensity, intref, atol=1e-5, rtol=0.02)

def test_antiferro_chain():
    unit_cell = UnitCell(3, 8, 8)
    #sites = generate_helical_structure(unit_cell, positions=[[0,0,0]], spins=[[0,1,0]],
    #                                   perpendicular=[0,0,1], propagation_vector=[0.5, 0, 0], names=["MCu1"])
    sites = generate_structure(unit_cell, positions=[[0,0,0]], spins=[[0, 1, 0]],
                               propagation_vectors=[0.5, 0, 0], names=["MCu1"])
    exchanges = generate_exchanges(sites=sites, max_distance=3.1, j=1)
    hamiltonian = Hamiltonian(sites, exchanges)
    energy, intensity = hamiltonian.energies_and_intensities(QVECS, use_rust=False, intensity_unit='spin')
    _test_ref_data('antiferro_chain', *omegasum(energy, intensity, is_series=False))

def test_antiferro_with_field():
    unit_cell = UnitCell(4,6,6)
    sites = generate_structure(unit_cell, positions=[[0,0,0], [0.5,0,0]], spins=[[0, 0, 1], [0, 0, -1]], names=['X', 'Y'])
    exchanges = generate_exchanges(sites=sites, bond=1, j=1, direction_filter=filter([1,0,0], symmetric=True))
    anisotropies = axis_anisotropies(sites, -0.1)
    hamiltonian = Hamiltonian(sites, exchanges, anisotropies)
    energy, intensity = hamiltonian.energies_and_intensities(QVECS, field=[0,0,7], use_rust=False)
    _test_ref_data('antiferro_chain_with_field', *omegasum(energy, intensity, is_series=False))

def test_kagome_antiferro():
    unit_cell = UnitCell(6, 6, 10, gamma=120)
    s = generate_structure(unit_cell, positions=[[0.5,0,0], [0,0.5,0], [0.5,0.5,0]], spins=[[1, 2, 0], [-2, -1, 0], [1, -1, 0]],
                           names=['X','Y','Z'], magnitudes=[1,1,1], spins_unit='lu')
    j1 = generate_exchanges(sites=s, bond=1, j=1)
    j2 = generate_exchanges(sites=s, bond=2, j=0.11)
    exchanges = j1 + j2
    hamiltonian = Hamiltonian(s, exchanges)
    energy, intensity = hamiltonian.energies_and_intensities(QVECS, use_rust=False, use_rotating=False)
    _test_ref_data('kagome_antiferro', *omegasum(energy, intensity, is_series=False))

def test_kagome_supercell():
    unit_cell = UnitCell(6, 6, 40, gamma=120)
    s = generate_helical_structure(unit_cell, positions=[[0.5,0,0], [0,0.5,0], [0.5,0.5,0]],
                                   spins=[[0, 1, 0], [0, 1, 0], [-1, -1, 0]], magnitudes=[1, 1, 1], names=['X', 'Y', 'Z'],
                                   spins_unit='lu', perpendicular=[0,0,1], propagation_vector=[-1./3., -1./3., 0])
    exchanges = generate_exchanges(sites=s, bond=1, j=1)
    hamiltonian = Hamiltonian(s, exchanges)
    energy, intensity = hamiltonian.energies_and_intensities(QVECS, use_rust=False, use_rotating=False)
    _test_ref_data('kagome_supercell', *omegasum(energy, intensity, is_series=False))

def test_tri_antiferro():
    unit_cell = UnitCell(3, 3, 4, gamma=120)
    sites = generate_helical_structure(unit_cell, positions=[[0,0,0]], spins=[[0, 1, 0]], magnitudes=[3. / 2], names=['X'],
                                       perpendicular=[0,0,1], propagation_vector=[1./3., 1./3., 0])
    exchanges = generate_exchanges(sites=sites, bond=1, j=1)
    anisotropies = axis_anisotropies(sites, 0.2)
    hamiltonian = Hamiltonian(sites, exchanges, anisotropies)
    energy, intensity = hamiltonian.energies_and_intensities(QVECS, use_rust=False, use_rotating=True)
    # Ignore last Q point as it is not positive definite and numerically different in Matlab and Python
    _test_ref_data('tri_antiferro', *omegasum(energy, intensity, is_series=False), ignoreQ=[5])

def test_square_antiferro():
    unit_cell = UnitCell(3, 3, 9)
    sites = generate_structure(unit_cell, positions=[[0,0,0]], spins=[[1,0,0]], names=['X'],
                               propagation_vectors=[[0.5, 0.5, 0]])
    # Note this model has spins in a-b plane but z-axis anisotropy so will give imaginary modes (needs non-herm solver)
    hamiltonian = Hamiltonian(sites, generate_exchanges(sites, bond=1, j=0.5), axis_anisotropies(sites, -0.1))
    energy, intensity = hamiltonian.energies_and_intensities(QVECS, use_rust=False, use_rotating=False)
    _test_ref_data('square_antiferro', *omegasum(energy, intensity, is_series=False))
