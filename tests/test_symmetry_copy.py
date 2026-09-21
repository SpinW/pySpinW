""" Tests for symmetry_copy methods """
import numpy as np
import pytest

from pyspinw.exchange import Exchange
from pyspinw.anisotropy import Anisotropy

from symmetry.exchange_test_cases import all_cases
from symmetry.anisotropy_test_cases import anisotropy_test_system

@pytest.mark.parametrize("hamiltonian", all_cases)
def test_symmetry_copy_against_symmetry_fill_exchanges(hamiltonian):
    """ Makes sure symmetry copy and symmetry fill are consistent for exchanges

    Note: this means that if symmetry_fill is wrong, then a passing test means symmetry_copy is also wrong
    """

    for exchange in hamiltonian.exchanges:
        others: list[Exchange] = exchange.symmetry_fill(hamiltonian.structure)

        for to_replicate in others:
            site_1 = to_replicate.site_1
            site_2 = to_replicate.site_2
            offset = to_replicate.cell_offset

            replicated = exchange.symmetry_copy(hamiltonian.structure, site_1, site_2, offset)

            assert np.allclose(to_replicate.exchange_matrix, replicated.exchange_matrix)


def test_symmetry_copy_against_symmetry_fill_anisotropies():
    """ Makes sure symmetry copy and symmetry fill are consistent for anisotropies

    Note: this means that if symmetry_fill is wrong, then a passing test means symmetry_copy is also wrong
    """

    hamiltonian = anisotropy_test_system()

    for anisotropy in hamiltonian.anisotropies:
        others: list[Anisotropy] = anisotropy.symmetry_fill(hamiltonian.structure)

        for to_replicate in others:
            site = to_replicate.site

            replicated = anisotropy.symmetry_copy(hamiltonian.structure, site)

            assert np.allclose(to_replicate.anisotropy_matrix, replicated.anisotropy_matrix)