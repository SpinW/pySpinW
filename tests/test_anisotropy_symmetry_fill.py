""" Tests for anisotropy filling in Hamiltonian.symmetry_fill"""
from collections import Counter

import numpy as np
import pytest

from symmetry.anisotropy_test_cases import anisotropy_test_system

hamiltonian = anisotropy_test_system()
filled = hamiltonian.symmetry_filled()

def test_anisotropy_symmetry_fill_fills():
    """ Check that symmetry_filled fills all sites in the test example """

    # Check that every anisotropy matches exactly one site
    counter = Counter([anisotropy.site.unique_id for anisotropy in filled.anisotropies])

    for site in filled.structure.site:
        assert counter[site.unique_id] == 1

