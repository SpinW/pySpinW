import numpy as np
import pytest

import pyspinw
from pyspinw import *
from pyspinw.polarisation import calculate_polarised_intensity
from pyspinw.components.components import Mab, Sperp

def make_ferromagnetic_chain():
    """ A spin-1 ferromagnetic chain along a """
    site = Site(0, 0, 0, 0, 0, 1, name="S")
    coupling = HeisenbergExchange(site, site, j=-1.0, cell_offset=(1, 0, 0), name="J")
    return Hamiltonian(Structure([site], UnitCell(3.0, 8.0, 8.0)), [coupling])

chain = make_ferromagnetic_chain()
simple_path = Path([[0, 0, 0], [1, 0, 0]], n_points_per_segment=11)
less_simple_path = Path([[-1, -1, 0], [0, 0, 0], [0, 0, 1]], n_points_per_segment=11)

test_systems = [(chain, simple_path),
                (chain, less_simple_path)]

xyz = [["xx", "xy", "xz"], ["yx", "yy", "yz"], ["zx", "zy", "zz"]]



@pytest.mark.parametrize("i", range(3))
@pytest.mark.parametrize("j", range(3))
@pytest.mark.parametrize("use_rust", [True, False])
@pytest.mark.parametrize("hamiltonian, path", test_systems)
def test_mab_consistent(use_rust, i, j, hamiltonian, path):
    """ components='all' gives the (Pab, Pa, Mab) tuple, with q along the last axis """

    component = "M" + xyz[i][j]

    q_vectors = path.q_points()

    # Get the Sab matrix and other stuff
    energies, intensities, sab, _, scaling, rlu_to_cart = chain._spinwave_calculation(
        q_vectors=q_vectors,
        use_rust=use_rust,
        use_rotating=False,
        intensity_unit='cell',
        save_sab=True,
    )

    # Data from the original method
    original = calculate_polarised_intensity(sab,
        Sperp = intensities,
        q_vectors=q_vectors,
        components=component,
        rlu_to_cart=rlu_to_cart).T



    # Test the big matrix
    m_matrix = Mab(sab, q_vectors, rlu_to_cart, None)

    test_1 = m_matrix[i,j,:,:]

    print("Test data")
    print(test_1)

    # Test the MXX type

    globals().update({
        "test_component": getattr(pyspinw.components.components, component)
    })


    # IDE wont see this variable because we added it via globals()
    test_2 = test_component(sab, q_vectors, rlu_to_cart, None)

    difference = original - test_1
    assert np.allclose(difference, 0), f"M[{i}][{j}] should agree with original"

    difference = test_1 - test_2
    assert np.allclose(difference, 0), f"M[{i}][{j}] should agree with M{xyz[i][j]}"

    difference = original - test_2
    assert np.allclose(difference, 0), f"{xyz[i][j]} components should match original"

@pytest.mark.parametrize("hamiltonian, path", test_systems)
def test_sperp_consistent(hamiltonian, path):
    q_vectors = path.q_points()

    # Get the Sab matrix and other stuff
    energies, intensities, sab, _, scaling, rlu_to_cart = chain._spinwave_calculation(
        q_vectors=q_vectors,
        use_rotating=False,
        intensity_unit='cell',
        save_sab=True,
    )

    print("intensities:", intensities)

    # IDE wont see this variable because we added it via globals()
    sperp = Sperp(sab, q_vectors, rlu_to_cart, None)



    assert np.allclose(sperp, np.array(intensities).T)



@pytest.mark.parametrize("direction", ["x", "y", "z"])
@pytest.mark.parametrize("hamiltonian, path", test_systems)
def test_half_polarisation_consistent(hamiltonian, path, direction):
    q_vectors = path.q_points()

    # Get the Sab matrix and other stuff
    energies, intensities, sab, _, scaling, rlu_to_cart = chain._spinwave_calculation(
        q_vectors=q_vectors,
        use_rotating=False,
        intensity_unit='cell',
        save_sab=True,
    )

    component = "P" + direction
    globals().update({
        "test_component": getattr(pyspinw.components.components, component)
    })

    original = calculate_polarised_intensity(sab,
        Sperp = intensities,
        q_vectors=q_vectors,
        components=component,
        rlu_to_cart=rlu_to_cart).T

    # IDE wont see this variable because we added it via globals()
    pa = test_component(sab, q_vectors, rlu_to_cart, None)

    print("original", original)
    print("new", pa)

    assert np.allclose(pa, original)



