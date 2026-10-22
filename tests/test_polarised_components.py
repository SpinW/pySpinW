""" Tests for the components argument of Hamiltonian.energies_and_intensities """

import numpy as np
import pytest

from pyspinw import HeisenbergExchange, Hamiltonian, Site, Path, Structure, UnitCell


@pytest.fixture
def ferromagnetic_chain():
    """ A spin-1 ferromagnetic chain along a """
    site = Site(0, 0, 0, 0, 0, 1, name="S")
    coupling = HeisenbergExchange(site, site, j=-1.0, cell_offset=(1, 0, 0), name="J")
    return Hamiltonian(Structure([site], UnitCell(3.0, 8.0, 8.0)), [coupling])


@pytest.fixture
def path():
    """ Eleven points from (0,0,0) to (1,0,0) """
    return Path([[0, 0, 0], [1, 0, 0]], n_points_per_segment=11)


@pytest.mark.parametrize("use_rust", [True, False])
def test_all_components_returns_full_polarisation_set(ferromagnetic_chain, path, use_rust):
    """ components='all' gives the (Pab, Pa, Mab) tuple, with q along the last axis """
    n_q, n_modes = 11, 2

    energies, intensities = ferromagnetic_chain.energies_and_intensities(path, use_rust=use_rust, components="all")

    assert np.shape(energies) == (n_q, n_modes)
    assert isinstance(intensities, tuple)
    pab, pa, mab = intensities
    assert pab.shape == (3, 3, n_modes, n_q)
    assert pa.shape == (3, n_modes, n_q)
    assert mab.shape == (3, 3, n_modes, n_q)

    # The tuple is the same calculation as a single component, kept complex
    _, pxx = ferromagnetic_chain.energies_and_intensities(path, use_rust=use_rust, components="Pxx")
    np.testing.assert_allclose(pab[0, 0].real.T, pxx)
    assert np.any(mab.imag != 0), "Mab of a precessing ferromagnet has imaginary off-diagonal terms"


def test_unparseable_components_raises_value_error(ferromagnetic_chain, path):
    """ A string that is not a polarisation component expression names the accepted forms """
    with pytest.raises(ValueError, match="Pab"):
        ferromagnetic_chain.energies_and_intensities(path, use_rust=False, components="Sxx")


@pytest.mark.parametrize("leading_sign, factor", [("-", -1.0), ("+", 1.0)])
def test_leading_sign_is_the_prefactor_of_the_first_term(ferromagnetic_chain, path, leading_sign, factor):
    """ A sign before the first term multiplies it, as a sign before any later term does """
    _, pxx = ferromagnetic_chain.energies_and_intensities(path, use_rust=False, components="Pxx")
    _, pyy = ferromagnetic_chain.energies_and_intensities(path, use_rust=False, components="Pyy")
    assert np.any(pxx != 0), "the check needs a nonzero first term"

    _, combined = ferromagnetic_chain.energies_and_intensities(path, use_rust=False,
                                                               components=f"{leading_sign}Pxx+Pyy")

    np.testing.assert_allclose(combined, factor * pxx + pyy)
