"""Check the energy-gain branch of a non-reciprocal model against SpinW.

By detailed balance S(Q, -w) = exp(-w/kT) S(-Q, w), so the energy-gain mode at q is the loss mode
at -q: the negative branch at q must be -E(-q), not -E(q). The reference values are from SpinW 4
(MATLAB) `spinwave` + `sw_neutron` on the same chains.
"""

import numpy as np
import pytest

from pyspinw import Hamiltonian, HeisenbergExchange, LatticeSite, Path, Structure, UnitCell
from pyspinw.anisotropy import AxisMagnitudeAnisotropy
from pyspinw.exchange import DMExchange

# The CI job that builds the Rust core runs only tests marked `rust`
USE_RUST = [pytest.param(True, marks=pytest.mark.rust), False]


def dm_chain(dm: float, hard_x: float = 0.0):
    """Ferromagnetic chain with a DM vector along the spins, optionally with a hard axis along x."""
    site = LatticeSite(0, 0, 0, 0, 0, 1, name="S")
    exchanges = [
        HeisenbergExchange(site, site, j=-1.0, cell_offset=(1, 0, 0), name="J"),
        DMExchange(site, site, 0.0, 0.0, dm, cell_offset=(1, 0, 0), name="D"),
    ]
    anisotropies = [AxisMagnitudeAnisotropy(site, -0.06, direction=[0, 0, 1])]
    if hard_x:
        anisotropies.append(AxisMagnitudeAnisotropy(site, hard_x, direction=[1, 0, 0]))
    return Hamiltonian(Structure([site], UnitCell(3.0, 8.0, 8.0)), exchanges, anisotropies)


def ends_of_path(hamiltonian: Hamiltonian, use_rust: bool):
    """Energies and intensities at h = +0.1 and h = -0.1, each sorted by energy."""
    # 21 points, so that use_rust=False has more q-points than chunks
    path = Path([[0.1, 0, 0], [-0.1, 0, 0]], n_points_per_segment=21, avoid_endpoints=False)
    energies, intensities = hamiltonian.energies_and_intensities(path, use_rust=use_rust)
    energies, intensities = np.real(np.asarray(energies)), np.real(np.asarray(intensities))
    order = np.argsort(energies, axis=1)
    energies, intensities = np.take_along_axis(energies, order, 1), np.take_along_axis(intensities, order, 1)
    return energies[[0, -1]], intensities[[0, -1]]


@pytest.mark.parametrize("use_rust", USE_RUST)
@pytest.mark.parametrize("hard_x, energy_plus, energy_minus, intensity", [
    # B = 0: the loss branch is right without the fix, only the gain branch moves
    (0.0, [-0.795859, 0.208073], [-0.208073, 0.795859], 0.5),
    # B != 0: the lower right block also couples into the loss branch and the intensities
    (0.1, [-0.887494, 0.299709], [-0.299709, 0.887494], 0.591277),
])
def test_gain_branch_is_minus_loss_branch_at_minus_q(use_rust, hard_x, energy_plus, energy_minus, intensity):
    """The negative branch at q is -E(-q), as SpinW gives."""
    energies, intensities = ends_of_path(dm_chain(0.25, hard_x), use_rust)

    np.testing.assert_allclose(energies[0], energy_plus, atol=2e-6)
    np.testing.assert_allclose(energies[1], energy_minus, atol=2e-6)
    np.testing.assert_allclose(intensities, intensity, atol=2e-6)


@pytest.mark.parametrize("use_rust", USE_RUST)
def test_reciprocal_model_keeps_mirrored_branches(use_rust):
    """Without DM, E(q) = E(-q), so the gain branch at q stays the mirror of the loss branch."""
    energies, intensities = ends_of_path(dm_chain(0.0), use_rust)

    np.testing.assert_allclose(energies[0], -energies[0][::-1], atol=1e-9)
    np.testing.assert_allclose(energies[0], energies[1], atol=1e-9)
    np.testing.assert_allclose(intensities[0], intensities[1], atol=1e-9)
