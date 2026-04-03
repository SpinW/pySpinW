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

def _test_ref_data(test_name, energy, intensity):
    eref, intref = omegasum(REFDAT[test_name][0][0].T, REFDAT[test_name][0][1].T, is_series=False)
    # Remove very large intensities as this diverges near zone centre
    intensity[np.where(intensity > INTHIGH)] = np.nan 
    intref[np.where(intref > INTHIGH)] = np.nan 
    np.testing.assert_allclose(energy, eref, atol=1e-5, rtol=0)
    np.testing.assert_allclose(intensity, intref, atol=1e-5, rtol=0.02)

def test_antiferro_chain():
    unit_cell = UnitCell(3, 8, 8)
    #sites = generate_helical_structure(unit_cell, positions=[[0,0,0]], moments=[[0,1,0]],
    #                                   perpendicular=[0,0,1], propagation_vector=[0.5, 0, 0], names=["MCu1"])
    sites = generate_structure(unit_cell, positions=[[0,0,0]], moments=[[0,1,0]],
                               propagation_vectors=[0.5, 0, 0], names=["MCu1"])
    exchanges = generate_exchanges(sites=sites, max_distance=3.1, j=1)
    hamiltonian = Hamiltonian(sites, exchanges)
    energy, intensity = hamiltonian.energies_and_intensities(QVECS, use_rust=False, intensity_unit='spin')
    _test_ref_data('antiferro_chain', *omegasum(energy, intensity, is_series=False))

def test_antiferro_with_field():
    unit_cell = UnitCell(4,6,6)
    sites = generate_structure(unit_cell, positions=[[0,0,0], [0.5,0,0]], moments=[[0,0,1], [0,0,-1]], names=['X', 'Y'])
    exchanges = generate_exchanges(sites=sites, bond=1, j=1, direction_filter=filter([1,0,0], symmetric=True))
    anisotropies = axis_anisotropies(sites, -0.1)
    hamiltonian = Hamiltonian(sites, exchanges, anisotropies)
    energy, intensity = hamiltonian.energies_and_intensities(QVECS, field=[0,0,7], use_rust=False)
    _test_ref_data('antiferro_chain_with_field', *omegasum(energy, intensity, is_series=False))

def test_kagome_antiferro():
    unit_cell = UnitCell(6, 6, 10, gamma=120)
    s = generate_structure(unit_cell, positions=[[0.5,0,0], [0,0.5,0], [0.5,0.5,0]], moments=[[1,2,0], [-2,-1,0], [1,-1,0]],
                           names=['X','Y','Z'], magnitudes=[1,1,1], moments_unit='lu')
    j1 = generate_exchanges(sites=s, bond=1, j=1)
    j2 = generate_exchanges(sites=s, bond=2, j=0.11)
    exchanges = j1 + j2
    hamiltonian = Hamiltonian(s, exchanges)
    energy, intensity = hamiltonian.energies_and_intensities(QVECS, use_rust=False, use_rotating=False)
    _test_ref_data('kagome_antiferro', *omegasum(energy, intensity, is_series=False))

def test_kagome_supercell():
    unit_cell = UnitCell(6, 6, 40, gamma=120)
    s = generate_helical_structure(unit_cell, positions=[[0.5,0,0], [0,0.5,0], [0.5,0.5,0]],
                                   moments=[[0,1,0], [0,1,0], [-1,-1,0]], magnitudes=[1,1,1], names=['X', 'Y', 'Z'],
                                   moments_unit='lu', perpendicular=[0,0,1], propagation_vector=[-1./3., -1./3., 0])
    exchanges = generate_exchanges(sites=s, bond=1, j=1)
    hamiltonian = Hamiltonian(s, exchanges)
    energy, intensity = hamiltonian.energies_and_intensities(QVECS, use_rust=False, use_rotating=False)
    _test_ref_data('kagome_supercell', *omegasum(energy, intensity, is_series=False))

def test_tri_antiferro():
    unit_cell = UnitCell(3, 3, 4, gamma=120)
    sites = generate_helical_structure(unit_cell, positions=[[0,0,0]], moments=[[0,1,0]], magnitudes=[3./2], names=['X'],
                                       perpendicular=[0,0,1], propagation_vector=[1./3., 1./3., 0])
    exchanges = generate_exchanges(sites=sites, bond=1, j=1)
    anisotropies = axis_anisotropies(sites, 0.2)
    hamiltonian = Hamiltonian(sites, exchanges, anisotropies)
    energy, intensity = hamiltonian.energies_and_intensities(QVECS, use_rust=False, use_rotating=True)
    _test_ref_data('tri_antiferro', *omegasum(energy, intensity, is_series=False))
