""" Powder averages are over a sphere of fixed |q| in inverse angstrom (issue #210) """

import numpy as np
import pytest

from pyspinw.exchange import HeisenbergExchange
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.path import Path1D
from pyspinw.sample import Powder
from pyspinw.site import LatticeSite
from pyspinw.structure import Structure
from pyspinw.symmetry.unitcell import UnitCell


def ferromagnetic_chain(a: float) -> Hamiltonian:
    """ Chain along x with J = -1, S = 1, so E(h) = 2 (1 - cos 2 pi h) meV, band top 4 meV """
    site = LatticeSite(0, 0, 0, 0, 0, 1, name="S")
    exchange = HeisenbergExchange(site, site, j=-1.0, cell_offset=(1, 0, 0))
    return Hamiltonian(Structure([site], UnitCell(a, 8.0, 8.0)), [exchange])


def spectrum_at(a: float, q: float, max_energy: float, n_energy_bins: int) -> np.ndarray:
    """ Powder spectrum of the chain with lattice constant a at a single |q| in inverse angstrom """
    path = Path1D(q, q, n_q=1, avoid_endpoints=False)
    return Powder(ferromagnetic_chain(a)).spectrum(path, n_samples=400, min_energy=0.0, max_energy=max_energy,
                                                   n_energy_bins=n_energy_bins)[-1][0]


@pytest.mark.parametrize("a", [2.0, 3.0])
def test_no_powder_intensity_above_the_highest_energy_at_this_q(a):
    """ Check that the powder spectrum is not found at energies that are too high.

    Along the chain, the largest component of a vector of length |q| in inverse angstrom is
    h = |q| a / 2 pi in lattice units, so the highest magnon energy in the powder average is E(h).
    There should be intensity below that and none above it. Averaging over a sphere of radius |q|
    in lattice units instead would reach h = |q| and give intensity up to E(|q|), which is higher
    for these cells.
    """
    q = 0.6
    h_max = q * a / (2 * np.pi)
    e_max = 2 * (1 - np.cos(2 * np.pi * h_max))

    # Bins of 0.05 meV, up to E(|q|) so that the range covers where the old behaviour put intensity
    e_top = 2 * (1 - np.cos(2 * np.pi * q))
    n_bins = int(np.ceil(e_top / 0.05))
    spectrum = spectrum_at(a, q, max_energy=n_bins * 0.05, n_energy_bins=n_bins)
    bin_lower_edges = np.arange(n_bins) * 0.05

    assert np.sum(spectrum[bin_lower_edges <= e_max]) > 0
    assert np.all(spectrum[bin_lower_edges > e_max + 0.05] == 0)


def test_spectrum_depends_on_q_times_lattice_constant():
    """ Check that the powder spectrum depends on |q| and the lattice constant only through |q| a.

    Scaling the lattice constant up by some factor and |q| down by the same factor leaves
    h = |q| a / 2 pi, and so the spectrum, unchanged.
    """
    q = 0.3
    spectrum = spectrum_at(2 * np.pi, q, max_energy=4.0, n_energy_bins=40)
    reference = spectrum_at(1.0, q * 2 * np.pi, max_energy=4.0, n_energy_bins=40)

    assert spectrum == pytest.approx(reference, rel=1e-9, abs=1e-12)
