""" Tests that the public calculation methods take an array of q-vectors as well as a Path """

import numpy as np
import pytest

from pyspinw import HeisenbergExchange, Hamiltonian, LatticeSite, Path, Structure, UnitCell


@pytest.fixture
def ferromagnetic_chain():
    """ A spin-1 ferromagnetic chain along a """
    site = LatticeSite(0, 0, 0, 0, 0, 1, name="S")
    coupling = HeisenbergExchange(site, site, j=-1.0, cell_offset=(1, 0, 0), name="J")
    return Hamiltonian(Structure([site], UnitCell(3.0, 8.0, 8.0)), [coupling])


@pytest.fixture
def path():
    """ A path through three points that are not along one line, with its endpoints included """
    return Path([[0.1, 0, 0], [0.3, 0.2, 0], [0.25, 0, 0.5]], n_points_per_segment=5, avoid_endpoints=False)


@pytest.mark.parametrize("use_rust", [True, False])
def test_energies_and_intensities_takes_an_array(ferromagnetic_chain, path, use_rust):
    """ An n-by-3 array gives the same result as a Path through the same q-points """
    q_vectors = path.q_points().tolist()

    from_array = ferromagnetic_chain.energies_and_intensities(q_vectors, use_rust=use_rust)

    from_path = ferromagnetic_chain.energies_and_intensities(path, use_rust=use_rust)
    for array_result, path_result in zip(from_array, from_path):
        np.testing.assert_array_equal(array_result, path_result)


@pytest.mark.parametrize("use_rust", [True, False])
def test_spinwave_calculation_takes_an_array(ferromagnetic_chain, path, use_rust):
    """ An n-by-3 array gives the same result as a Path through the same q-points """
    q_vectors = path.q_points()

    energies, intensities, _, _ = ferromagnetic_chain.spinwave_calculation(q_vectors, use_rust=use_rust)

    path_energies, path_intensities, _, _ = ferromagnetic_chain.spinwave_calculation(path, use_rust=use_rust)
    np.testing.assert_array_equal(energies, path_energies)
    np.testing.assert_array_equal(intensities, path_intensities)


@pytest.mark.parametrize("method", ["energies_and_intensities", "spinwave_calculation"])
@pytest.mark.parametrize("q_vectors", [[0.1, 0, 0], [[0.1, 0], [0.2, 0]]])
def test_wrong_array_shape_raises(ferromagnetic_chain, method, q_vectors):
    """ An array that is not n-by-3 is refused with its shape named """
    with pytest.raises(ValueError, match="n-by-3"):
        getattr(ferromagnetic_chain, method)(q_vectors, use_rust=False)
