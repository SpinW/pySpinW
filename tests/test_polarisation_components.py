import numpy as np
import pytest

from pyspinw import *
from pyspinw.polarisation import calculate_polarised_intensity
from pyspinw.components.components import Mab

def make_ferromagnetic_chain():
    """ A spin-1 ferromagnetic chain along a """
    site = Site(0, 0, 0, 0, 0, 1, name="S")
    coupling = HeisenbergExchange(site, site, j=-1.0, cell_offset=(1, 0, 0), name="J")
    return Hamiltonian(Structure([site], UnitCell(3.0, 8.0, 8.0)), [coupling])

chain = make_ferromagnetic_chain()

path = Path([[0, 0, 0], [1, 0, 0]], n_points_per_segment=11)

xyz = [["xx", "xy", "xz"], ["yx", "yy", "yz"], ["zx", "zy", "zz"]]

@pytest.mark.parametrize("i", range(3))
@pytest.mark.parametrize("j", range(3))
@pytest.mark.parametrize("use_rust", [True, False])
def test_mab_consistent(use_rust, i, j):
    """ components='all' gives the (Pab, Pa, Mab) tuple, with q along the last axis """



    q_vectors = path.q_points()

    # Get the Sab matrix and other stuff
    energies, intensities, sab, _, scaling, rlu_to_cart = chain._spinwave_calculation(
        q_vectors=q_vectors,
        use_rust=use_rust,
        use_rotating=False,
        intensity_unit='cell',
        save_sab=True,
    )


    m_matrix = Mab(sab, q_vectors, rlu_to_cart, None)

    for i in range(3):
        for j in range(3):

            original = calculate_polarised_intensity(sab,
                Sperp = intensities,
                q_vectors=q_vectors,
                components="M"+xyz[j][i],
                rlu_to_cart=rlu_to_cart)

            test = m_matrix[i,j,:,:]

            difference = original - test.T

            assert np.allclose(difference, 0), f"{xyz[i][j]} components should match"

