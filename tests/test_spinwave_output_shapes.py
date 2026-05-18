import pytest

from pyspinw import *


@pytest.mark.parametrize("use_rust", [True, False])
def test_shape_sizes_non_rotating(use_rust: bool):
    unit_cell = UnitCell(6, 6, 5, gamma=120)

    x = LatticeSite(0.5, 0, 0, 0, 1, 0, name="X")
    y = LatticeSite(0, 0.5, 0, 0, 1, 0, name="Y")
    z = LatticeSite(0.5, 0.5, 0, 0, 1, 0, name="Z")

    sites = [x, y, z]

    s = Structure(sites, unit_cell=unit_cell, supercell=TiledSupercell(scaling=(1, 1, 1)))

    exchanges = generate_exchanges(sites=[x, y, z],
                                   unit_cell=unit_cell,
                                   max_distance=4.,
                                   exchange_type=HeisenbergExchange,
                                   j=-1)

    hamiltonian = Hamiltonian(s, exchanges)

    path = Path([[-0.5, 0, 0], [0, 0, 0], [0.5, 0.5, 0]])

    q = path.q_points()

    n_points = q.shape[0]
    n_sites = len(hamiltonian.expanded().structure.sites)


    energies, intensities, sab, wavefunctions, _, _ = hamiltonian._spinwave_calculation(q, field=None, use_rotating=False,
                                      use_rust=use_rust,
                                      save_sab=True, save_wavefunctions=True)

    # Check energy
    assert len(energies) == n_points, "There should be the same number of energy entries as q points"
    assert all([len(levels) == 2*n_sites for levels in energies]), "Each energy result should have 2n_sites entries"

    # Check intensity
    assert len(intensities) == n_points, "There should be the same number of intensity entries as q points"
    assert all([len(levels) == 2 * n_sites for levels in intensities]), (
        "Each intensity result should have 2n_sites entries")

    # Check Sab
    assert sab.shape == (3, 3, 2*n_sites, n_points)

    # Check wavefunctions
    assert len(wavefunctions) == n_points, "There should be the same number of wavefunction entries as q points"
    assert all([wavefunction.shape == (2*n_sites, 2*n_sites) for wavefunction in wavefunctions])

