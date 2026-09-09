""" Tests for anisotropy filling in Hamiltonian.symmetry_fill"""
from collections import Counter

import numpy as np
import pytest

from pyspinw import Anisotropy
from symmetry.anisotropy_test_cases import anisotropy_test_system

hamiltonian = anisotropy_test_system()
filled = hamiltonian.symmetry_filled()

# Mapping from site to anisotropy, should be unique (first test tests this)
site_to_anisotropy = {anisotropy.site.unique_id: anisotropy for anisotropy in filled.anisotropies}

def test_anisotropy_symmetry_fill_fills():
    """ Check that symmetry_filled fills all sites in the test example """

    # Check that every anisotropy matches exactly one site
    counter = Counter([anisotropy.site.unique_id for anisotropy in filled.anisotropies])

    for site in filled.structure.sites:
        assert counter[site.unique_id] == 1



@pytest.mark.parametrize("anisotropy", filled.anisotropies, ids=lambda aniso: aniso.site.name)
def test_anisotropy_symmetry_fill_transforms(anisotropy: Anisotropy):
    """ Check that there are appropriately transformed anisotropies"""
    site = anisotropy.site

    for other_site, operations in filled.structure.symmetry_related(site):
        other_anisotropy = site_to_anisotropy[other_site.unique_id]
        for operation in operations:
            transform = operation.point_operation_in_cartesian(filled.structure.unit_cell)

            expected_matrix = transform @ anisotropy.anisotropy_matrix @ transform.T

            assert np.allclose(other_anisotropy.anisotropy_matrix, expected_matrix)




