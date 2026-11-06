"""Tests for the ``NeighbourSearch`` class in ``pyspinw.nearest_neighbours``.

This test module verifies that ``NeighbourSearch``:
1. Always terminates (guarded against infinite loops).
2. Raises appropriate exceptions on invalid inputs (empty ``search_sites``,
   negative ``n``) and handles boundary cases such as ``n = 0``.
3. Returns the exact ``0``-th through ``n``-th nearest-neighbour shells in
   periodic crystals (for ``n`` up to 10) both when ``start_site`` belongs to
   ``search_sites`` and when ``start_site`` does not belong to ``search_sites``.
4. Works across 1D (one unit-cell axis much shorter than the others),
   2D (one unit-cell axis much longer than the others), and 3D geometries.
5. Handles diverse unit-cell side lengths ``(a, b, c)`` and angles
   ``(alpha, beta, gamma)`` across all seven crystal systems, including
   orthogonal, hexagonal, rhombohedral, monoclinic, and acute/obtuse triclinic
   unit cells.
6. Works for both low-symmetry and high-symmetry crystals constructed via
   ``Structure`` and ``SpaceGroup``, with both small (1 site) and large
   (up to 48 sites) numbers of sites per unit cell.
"""

import math
import sys
import types
from collections import defaultdict
from unittest.mock import MagicMock

import numpy as np
import pytest

# Stub out the unused ``from plotly.express import histogram`` auto-import in
# ``pyspinw/nearest_neighbours.py`` if ``plotly`` is not installed in the test
# environment.
if "plotly" not in sys.modules:
    _plotly_mod = types.ModuleType("plotly")
    _plotly_express_mod = types.ModuleType("plotly.express")
    _plotly_express_mod.histogram = None
    _plotly_mod.express = _plotly_express_mod
    sys.modules["plotly"] = _plotly_mod
    sys.modules["plotly.express"] = _plotly_express_mod

from pyspinw.interface import get_spacegroup
from pyspinw.nearest_neighbours import NeighbourSearch
from pyspinw.site import Site
from pyspinw.structure import Structure
from pyspinw.symmetry.unitcell import UnitCell
from pyspinw.tolerances import tolerances
from pyspinw.util import IncrementalPointHistogram


MAX_OUTER_LOOP_ITERATIONS = 100

# ---------------------------------------------------------------------------
# Reference helper functions for verifying nearest-neighbour shells
# ---------------------------------------------------------------------------


def pair_distance(
    start_site: Site,
    target_site: Site,
    offset: tuple[int, int, int],
    unit_cell: UnitCell,
) -> float:
    """Compute the Euclidean distance between ``start_site`` and ``target_site + offset``.

    Parameters
    ----------
    start_site : Site
        Central reference site in cell ``(0, 0, 0)``.
    target_site : Site
        Neighbour site within the unit cell.
    offset : tuple of (int, int, int)
        Unit-cell lattice offset ``(di, dj, dk)`` applied to ``target_site``.
    unit_cell : UnitCell
        Crystal unit cell used to convert fractional coordinates to Cartesian.

    Returns
    -------
    float
        Cartesian distance in Angstroms.
    """
    delta_lu = np.asarray(target_site.ijk, dtype=float) + np.asarray(offset, dtype=float) - np.asarray(start_site.ijk, dtype=float)
    delta_cart = unit_cell.lattice_units_to_cartesian(delta_lu.reshape(1, 3))[0]
    return float(np.linalg.norm(delta_cart))


def brute_force_neighbour_shells(
    start_site: Site,
    search_sites: list[Site],
    unit_cell: UnitCell,
    n: int,
    shell_tolerance: float = tolerances.SAME_SITE_ABS_TOL,
) -> list[set[tuple[int, tuple[int, int, int]]]]:
    """Compute exact reference nearest-neighbour shells via supercell enumeration.

    Enumerates all unit-cell offsets ``(di, dj, dk)`` guaranteed to enclose a
    sphere containing at least ``n + 1`` shells around ``start_site``, groups
    candidates by squared Cartesian distance using ``IncrementalPointHistogram``
    with ``shell_tolerance``, and returns the ``0``-th through ``n``-th shells.

    Parameters
    ----------
    start_site : Site
        Site for which neighbours are searched.
    search_sites : list of Site
        Candidate neighbour sites in the unit cell.
    unit_cell : UnitCell
        Crystal unit cell.
    n : int
        Maximum nearest-neighbour shell order.
    shell_tolerance : float, optional
        Squared-distance grouping tolerance.

    Returns
    -------
    list of set of (int, tuple of (int, int, int))
        List of length ``n + 1`` where element ``k`` is the set of
        ``(site.unique_id, (di, dj, dk))`` pairs belonging to the ``k``-th
        nearest-neighbour shell.
    """
    start_in_search = start_site.unique_id in {s.unique_id for s in search_sites}
    if n == 0:
        if start_in_search:
            return [{(start_site.unique_id, (0, 0, 0))}]
        return [set()]

    # Determine a conservative physical radius guaranteed to contain >= n+1 shells.
    # Along the shortest lattice vector of length L_min, n+2 unit-cell steps
    # always produce at least n+2 distinct distances.
    basis = unit_cell.lattice_units_to_cartesian(np.eye(3, dtype=float))
    norms = np.linalg.norm(basis, axis=1)
    max_radius = float((n + 2) * np.min(norms) + np.max(norms))

    # Convert max_radius into required integer bounds (N_i, N_j, N_k) along each
    # reciprocal/perpendicular lattice direction using the metric tensor inverse.
    inv_basis = np.linalg.inv(basis)
    max_offsets = np.ceil(max_radius * np.linalg.norm(inv_basis, axis=0)).astype(int) + 2
    ni, nj, nk = (min(int(m), n + 4) for m in max_offsets)

    grid_i, grid_j, grid_k = np.meshgrid(
        np.arange(-ni, ni + 1),
        np.arange(-nj, nj + 1),
        np.arange(-nk, nk + 1),
        indexing="ij",
    )
    offsets = np.stack([grid_i.ravel(), grid_j.ravel(), grid_k.ravel()], axis=1)
    cart_offsets = unit_cell.lattice_units_to_cartesian(offsets.astype(float))

    start_pos = unit_cell.lattice_units_to_cartesian(np.asarray(start_site.ijk, dtype=float).reshape(1, 3))
    search_pos = unit_cell.lattice_units_to_cartesian(
        np.array([s.ijk for s in search_sites], dtype=float)
    )
    rel_pos = search_pos - start_pos  # (n_sites, 3)

    # Vectorised distance calculation across all candidate supercell offsets:
    # diffs has shape (n_offsets, n_sites, 3)
    diffs = rel_pos[np.newaxis, :, :] + cart_offsets[:, np.newaxis, :]
    dists_sq = np.sum(diffs ** 2, axis=2)

    # Sort all candidate pairs by squared distance and only populate the
    # histogram up to shell n + 1 so the reference check is O(N log N) and fast.
    flat_order = np.argsort(dists_sq, axis=None)
    n_sites = len(search_sites)
    target_groups = (n + 1) if start_in_search else n

    hist = IncrementalPointHistogram(tolerance=shell_tolerance)
    for flat_idx in flat_order:
        offset_idx, site_idx = divmod(int(flat_idx), n_sites)
        d_sq = float(dists_sq[offset_idx, site_idx])
        ox, oy, oz = offsets[offset_idx]
        entry = (search_sites[site_idx].unique_id, (int(ox), int(oy), int(oz)))
        hist.add(d_sq, entry)
        if hist.n_groups() > target_groups:
            break

    groups = hist.groups()
    if start_in_search:
        selected = groups[: n + 1]
    else:
        selected = [[]] + groups[:n]

    return [set(group) for group in selected]


def assert_shells_match_reference(
    shells: list[list[tuple[Site, tuple[int, int, int]]]],
    start_site: Site,
    search_sites: list[Site],
    unit_cell: UnitCell,
    n: int,
    shell_tolerance: float = tolerances.SAME_SITE_ABS_TOL,
) -> list[float]:
    """Assert that ``shells`` returned by ``NeighbourSearch.neighbour_sites()`` are valid and exact.

    Checks:
    1. Output has length ``n + 1``.
    2. Shell 0 contains ``[(start_site, (0, 0, 0))]`` if ``start_site`` is in
       ``search_sites``, or ``[]`` otherwise.
    3. Every shell ``k >= 1`` is non-empty, contains no duplicate
       ``(site.unique_id, offset)`` entries, and has a single well-defined
       Cartesian distance ``d_k`` shared by all members of the shell.
    4. Shell distances are strictly increasing: ``0 < d_1 < d_2 < ... < d_n``.
    5. Each shell ``k`` matches the exact reference shell from
       ``brute_force_neighbour_shells``.

    Returns
    -------
    list of float
        The characteristic Cartesian distance of each shell ``0..n`` (with
        ``NaN`` for an empty 0th shell when ``start_site`` is not in
        ``search_sites``).
    """
    assert len(shells) == n + 1, f"Expected {n + 1} shells (0..{n}), got {len(shells)}"

    start_in_search = start_site.unique_id in {s.unique_id for s in search_sites}
    if start_in_search:
        assert len(shells[0]) == 1, "0th shell should contain only start_site at (0, 0, 0)"
        s0, off0 = shells[0][0]
        assert s0.unique_id == start_site.unique_id
        assert tuple(off0) == (0, 0, 0)
    else:
        assert shells[0] == [], f"0th shell should be empty when start_site is not in search_sites, got {shells[0]}"

    shell_distances: list[float] = [0.0 if start_in_search else math.nan]
    for k in range(1, n + 1):
        shell = shells[k]
        assert len(shell) > 0, f"Shell {k} should not be empty"

        seen_pairs = set()
        dists_sq = []
        for site, offset in shell:
            pair_key = (site.unique_id, (int(offset[0]), int(offset[1]), int(offset[2])))
            assert pair_key not in seen_pairs, f"Duplicate site/offset {pair_key} in shell {k}"
            seen_pairs.add(pair_key)

            d = pair_distance(start_site, site, offset, unit_cell)
            dists_sq.append(d ** 2)

        assert max(dists_sq) - min(dists_sq) <= 2 * shell_tolerance, (
            f"Points in shell {k} do not share the same squared distance: "
            f"min={min(dists_sq)}, max={max(dists_sq)}"
        )
        shell_distances.append(math.sqrt(float(np.mean(dists_sq))))

    for k in range(1, n):
        assert shell_distances[k] < shell_distances[k + 1], (
            f"Shell distances are not strictly increasing: "
            f"d_{k}={shell_distances[k]} >= d_{k+1}={shell_distances[k+1]}"
        )

    expected_shells = brute_force_neighbour_shells(
        start_site, search_sites, unit_cell, n=n, shell_tolerance=shell_tolerance
    )
    for k in range(n + 1):
        actual_set = {
            (site.unique_id, (int(offset[0]), int(offset[1]), int(offset[2])))
            for site, offset in shells[k]
        }
        assert actual_set == expected_shells[k], (
            f"Mismatch in shell {k} (expected {len(expected_shells[k])} sites, "
            f"got {len(actual_set)} sites).\n"
            f"Missing: {expected_shells[k] - actual_set}\n"
            f"Extra:   {actual_set - expected_shells[k]}"
        )

    return shell_distances


# ===========================================================================
# 1. Input Validation and Edge / Boundary Cases
# ===========================================================================


def test_empty_search_sites_raises_value_error():
    """Verify that passing an empty ``search_sites`` list raises ``ValueError``.

    System description:
        A simple cubic unit cell ``a = b = c = 1`` with a valid ``start_site``
        at the origin, but ``search_sites = []``. Because no periodic lattice
        of target sites exists, ``NeighbourSearch`` must reject the input.
    """
    unit_cell = UnitCell(1.0, 1.0, 1.0)
    start_site = Site(0.0, 0.0, 0.0)

    with pytest.raises(ValueError, match="at least one search site"):
        NeighbourSearch(start_site, [], unit_cell, n=1, max_iters=MAX_OUTER_LOOP_ITERATIONS)


@pytest.mark.parametrize("negative_n", [-1, -2, -10])
def test_negative_n_raises_value_error(negative_n: int):
    """Verify that negative shell orders ``n < 0`` raise ``ValueError``.

    System description:
        A simple cubic unit cell ``a = b = c = 1`` with one site at ``(0, 0, 0)``,
        queried with negative integer values of ``n``.
    """
    unit_cell = UnitCell(1.0, 1.0, 1.0)
    site = Site(0.0, 0.0, 0.0)

    with pytest.raises(ValueError, match="positive"):
        NeighbourSearch(site, [site], unit_cell, n=negative_n, max_iters=MAX_OUTER_LOOP_ITERATIONS)


def test_n_zero_when_start_site_in_search_set():
    """Verify ``n = 0`` when ``start_site`` is included in ``search_sites``.

    System description:
        A simple cubic crystal (``a = b = c = 1``) with a single site at
        ``(0.5, 0.5, 0.5)``. Asking for ``n = 0`` should return a single-element
        outer list containing only the 0th shell ``[(start_site, (0, 0, 0))]``.
    """
    unit_cell = UnitCell(1.0, 1.0, 1.0)
    site = Site(0.5, 0.5, 0.5)

    search = NeighbourSearch(site, [site], unit_cell, n=0, max_iters=MAX_OUTER_LOOP_ITERATIONS)
    shells = search.neighbour_sites()

    assert len(shells) == 1
    assert shells[0] == [(site, (0, 0, 0))]


def test_n_zero_when_start_site_not_in_search_set():
    """Verify ``n = 0`` when ``start_site`` is NOT in ``search_sites``.

    System description:
        A simple cubic crystal (``a = b = c = 1``) where ``start_site`` is at
        ``(0, 0, 0)`` and ``search_sites`` contains a different site at
        ``(0.5, 0.5, 0.5)``. Because ``start_site`` is not in ``search_sites``,
        the 0th shell must be empty ``[]`` and ``neighbour_sites()`` must have
        length 1 (``[[]]``).
    """
    unit_cell = UnitCell(1.0, 1.0, 1.0)
    start_site = Site(0.0, 0.0, 0.0)
    other_site = Site(0.5, 0.5, 0.5)

    search = NeighbourSearch(start_site, [other_site], unit_cell, n=0, max_iters=MAX_OUTER_LOOP_ITERATIONS)
    shells = search.neighbour_sites()

    assert shells == [[]]


# ===========================================================================
# 2. 1D Periodic Crystal Systems (one unit-cell axis much shorter than others)
# ===========================================================================


@pytest.mark.parametrize(
    "chain_axis, cell_lengths, expected_step",
    [
        ("a", (0.2, 10.0, 10.0), (1, 0, 0)),
        ("b", (10.0, 0.2, 10.0), (0, 1, 0)),
        ("c", (10.0, 10.0, 0.2), (0, 0, 1)),
    ],
)
def test_1d_uniform_chain_start_in_search_set(
    chain_axis: str,
    cell_lengths: tuple[float, float, float],
    expected_step: tuple[int, int, int],
):
    """Verify 1st through 10th nearest neighbours on a uniform 1D chain along a, b, or c.

    System description:
        A 1D periodic chain simulated by setting one unit-cell axis length to
        ``0.2`` Angstrom and the two perpendicular axes to ``10.0`` Angstrom,
        with a single site at ``(0.5, 0.5, 0.5)``. Up to ``n = 10`` (maximum
        chain distance ``10 * 0.2 = 2.0 < 10.0``), all nearest neighbours lie
        strictly along the short axis at cell offsets ``±k * expected_step``
        with distance ``0.2 * k`` and multiplicity 2 per shell.
    """
    unit_cell = UnitCell(*cell_lengths)
    site = Site(0.5, 0.5, 0.5, name=f"chain_{chain_axis}")
    n = 10

    shells = NeighbourSearch(site, [site], unit_cell, n=n,
                             max_iters=MAX_OUTER_LOOP_ITERATIONS).neighbour_sites()

    distances = assert_shells_match_reference(shells, site, [site], unit_cell, n=n)

    sx, sy, sz = expected_step
    for k in range(1, n + 1):
        assert math.isclose(distances[k], 0.2 * k, rel_tol=1e-9, abs_tol=1e-9)
        offsets_k = {offset for _, offset in shells[k]}
        assert offsets_k == {(k * sx, k * sy, k * sz), (-k * sx, -k * sy, -k * sz)}


def test_1d_chain_with_transverse_crossover():
    """Verify the 1D chain example from ``nearest_neighbours.py`` where transverse shells appear.

    System description:
        Matches the ``__main__`` example: ``UnitCell(1, 0.2, 1)`` with one site
        at ``(0.5, 0.5, 0.5)`` and ``n = 10``.
        - Shells 1..4 are pure 1D neighbours along ``b`` at distances
          ``0.2, 0.4, 0.6, 0.8`` (2 sites each).
        - Shell 5 is at distance ``1.0``, where the 5th chain neighbours
          ``(0, ±5, 0)`` become degenerate with the 1st transverse neighbours
          ``(±1, 0, 0)`` and ``(0, 0, ±1)``, giving 6 sites.
        - Shells 6..10 include mixed ``(di, dj, dk)`` offsets.
    """
    unit_cell = UnitCell(1.0, 0.2, 1.0)
    site = Site(0.5, 0.5, 0.5)
    n = 10

    shells = NeighbourSearch(site, [site], unit_cell, n=n,
                             max_iters=MAX_OUTER_LOOP_ITERATIONS).neighbour_sites()
    distances = assert_shells_match_reference(shells, site, [site], unit_cell, n=n)

    # Shells 1..4: pure 1D along b
    for k in range(1, 5):
        assert len(shells[k]) == 2
        assert math.isclose(distances[k], 0.2 * k, abs_tol=1e-9)
        assert {off for _, off in shells[k]} == {(0, -k, 0), (0, k, 0)}

    # Shell 5: d = 1.0 has (0, ±5, 0), (±1, 0, 0), (0, 0, ±1) -> 6 neighbours
    assert len(shells[5]) == 6
    assert math.isclose(distances[5], 1.0, abs_tol=1e-9)
    assert {off for _, off in shells[5]} == {
        (0, -5, 0),
        (0, 5, 0),
        (-1, 0, 0),
        (1, 0, 0),
        (0, 0, -1),
        (0, 0, 1),
    }


def test_1d_chain_start_site_not_in_search_set():
    """Verify 1D chain shells when ``start_site`` is NOT in ``search_sites``.

    System description:
        A quasi-1D crystal ``UnitCell(10.0, 0.2, 10.0)`` along ``b``.
        ``start_site`` is at ``(0.5, 0.2, 0.5)`` and ``search_sites`` contains
        only ``target_site = Site(0.5, 0.0, 0.5)``.
        Because ``start_site`` is displaced by ``+0.2 * b`` from ``target_site``,
        the distances to ``target_site`` in cell ``(0, m, 0)`` are
        ``|m - 0.2| * 0.2``:
        - Shell 0: ``[]`` (``start_site`` not in ``search_sites``)
        - Shell 1: ``m = 0``  (distance ``0.2 * 0.2 = 0.04``)
        - Shell 2: ``m = +1`` (distance ``0.8 * 0.2 = 0.16``)
        - Shell 3: ``m = -1`` (distance ``1.2 * 0.2 = 0.24``)
        - Shell 4: ``m = +2`` (distance ``1.8 * 0.2 = 0.36``)
        - Shell 5: ``m = -2`` (distance ``2.2 * 0.2 = 0.44``)
        - Shell 6: ``m = +3`` (distance ``2.8 * 0.2 = 0.56``)
    """
    unit_cell = UnitCell(10.0, 0.2, 10.0)
    start_site = Site(0.5, 0.2, 0.5, name="interstitial")
    target_site = Site(0.5, 0.0, 0.5, name="chain_atom")
    n = 6

    shells = NeighbourSearch(start_site, [target_site], unit_cell, n=n,
                             max_iters=MAX_OUTER_LOOP_ITERATIONS).neighbour_sites()
    distances = assert_shells_match_reference(shells, start_site, [target_site], unit_cell, n=n)

    expected_m_sequence = [0, 1, -1, 2, -2, 3]
    for k, m in enumerate(expected_m_sequence, start=1):
        assert shells[k] == [(target_site, (0, m, 0))]
        assert math.isclose(distances[k], abs(m - 0.2) * 0.2, abs_tol=1e-9)


def test_1d_dimerized_two_site_chain():
    """Verify a 1D alternating-bond (dimerized) chain with two sites per unit cell.

    System description:
        A 1D chain along ``a`` with ``UnitCell(1.0, 20.0, 20.0)`` containing two
        sites per cell: ``s0`` at ``(0.0, 0.0, 0.0)`` and ``s1`` at
        ``(0.35, 0.0, 0.0)``.
        From ``s0``, the intra-cell distance to ``s1`` in ``(0, 0, 0)`` is
        ``0.35``, the inter-cell distance to ``s1`` in ``(-1, 0, 0)`` is
        ``0.65``, and the distance to ``s0`` in ``(±1, 0, 0)`` is ``1.00``.
        Tested up to ``n = 8`` both when searching ``[s0, s1]`` (``start_site``
        in ``search_sites``) and when searching only ``[s1]`` (``start_site``
        not in ``search_sites``).
    """
    unit_cell = UnitCell(1.0, 20.0, 20.0)
    s0 = Site(0.0, 0.0, 0.0, name="A")
    s1 = Site(0.35, 0.0, 0.0, name="B")
    n = 8

    # Case 1: start_site (s0) is in search_sites ([s0, s1])
    shells_all = NeighbourSearch(s0, [s0, s1], unit_cell, n=n,
                                 max_iters=MAX_OUTER_LOOP_ITERATIONS).neighbour_sites()
    dists_all = assert_shells_match_reference(shells_all, s0, [s0, s1], unit_cell, n=n)
    assert shells_all[1] == [(s1, (0, 0, 0))]
    assert math.isclose(dists_all[1], 0.35, abs_tol=1e-9)
    assert shells_all[2] == [(s1, (-1, 0, 0))]
    assert math.isclose(dists_all[2], 0.65, abs_tol=1e-9)
    assert {off for _, off in shells_all[3]} == {(-1, 0, 0), (1, 0, 0)}
    assert math.isclose(dists_all[3], 1.00, abs_tol=1e-9)

    # Case 2: start_site (s0) is NOT in search_sites ([s1] only)
    shells_b_only = NeighbourSearch(s0, [s1], unit_cell, n=n,
                                    max_iters=MAX_OUTER_LOOP_ITERATIONS).neighbour_sites()
    dists_b_only = assert_shells_match_reference(shells_b_only, s0, [s1], unit_cell, n=n)
    assert shells_b_only[0] == []
    assert shells_b_only[1] == [(s1, (0, 0, 0))]
    assert math.isclose(dists_b_only[1], 0.35, abs_tol=1e-9)
    assert shells_b_only[2] == [(s1, (-1, 0, 0))]
    assert math.isclose(dists_b_only[2], 0.65, abs_tol=1e-9)


# ===========================================================================
# 3. 2D Periodic Crystal Systems (one unit-cell axis much longer than others)
# ===========================================================================


def test_2d_square_lattice():
    """Verify coordination shells up to ``n = 8`` on a 2D square lattice.

    System description:
        A 2D square lattice in the ``ab``-plane simulated by setting
        ``a = b = 1.0`` and ``c = 10.0`` (``alpha = beta = gamma = 90``), with
        one site at ``(0, 0, 0)``.
        Expected 2D square-lattice shells ``k = 1..8``:
        - k=1: d^2 = 1  (±1, 0), (0, ±1)         -> 4 sites
        - k=2: d^2 = 2  (±1, ±1)                 -> 4 sites
        - k=3: d^2 = 4  (±2, 0), (0, ±2)         -> 4 sites
        - k=4: d^2 = 5  (±2, ±1), (±1, ±2)       -> 8 sites
        - k=5: d^2 = 8  (±2, ±2)                 -> 4 sites
        - k=6: d^2 = 9  (±3, 0), (0, ±3)         -> 4 sites
        - k=7: d^2 = 10 (±3, ±1), (±1, ±3)       -> 8 sites
        - k=8: d^2 = 13 (±3, ±2), (±2, ±3)       -> 8 sites
    """
    unit_cell = UnitCell(1.0, 1.0, 10.0)
    site = Site(0.0, 0.0, 0.0)
    n = 8

    shells = NeighbourSearch(site, [site], unit_cell, n=n,
                             max_iters=MAX_OUTER_LOOP_ITERATIONS).neighbour_sites()
    distances = assert_shells_match_reference(shells, site, [site], unit_cell, n=n)

    expected_sq_dists = [0, 1, 2, 4, 5, 8, 9, 10, 13]
    expected_sizes = [1, 4, 4, 4, 8, 4, 4, 8, 8]

    assert [len(s) for s in shells] == expected_sizes
    for k in range(1, n + 1):
        assert math.isclose(distances[k] ** 2, expected_sq_dists[k], abs_tol=1e-8)
        assert all(offset[2] == 0 for _, offset in shells[k])


def test_2d_triangular_lattice_non_orthogonal():
    """Verify coordination shells up to ``n = 6`` on a 2D triangular lattice (gamma = 120 deg).

    System description:
        A 2D hexagonal/triangular Bravais lattice in the ``ab``-plane simulated
        with ``UnitCell(1.0, 1.0, 10.0, gamma=120)`` and one site at ``(0, 0, 0)``.
        Expected 2D triangular-lattice shells ``k = 1..6``:
        - k=1: d = 1        -> 6 neighbours
        - k=2: d = sqrt(3)  -> 6 neighbours
        - k=3: d = 2        -> 6 neighbours
        - k=4: d = sqrt(7)  -> 12 neighbours
        - k=5: d = 3        -> 6 neighbours
        - k=6: d = sqrt(12) -> 6 neighbours
    """
    unit_cell = UnitCell(1.0, 1.0, 10.0, gamma=120.0)
    site = Site(0.0, 0.0, 0.0)
    n = 6

    shells = NeighbourSearch(site, [site], unit_cell, n=n,
                             max_iters=MAX_OUTER_LOOP_ITERATIONS).neighbour_sites()
    distances = assert_shells_match_reference(shells, site, [site], unit_cell, n=n)

    expected_sq_dists = [0, 1, 3, 4, 7, 9, 12]
    expected_sizes = [1, 6, 6, 6, 12, 6, 6]

    assert [len(s) for s in shells] == expected_sizes
    for k in range(1, n + 1):
        assert math.isclose(distances[k] ** 2, expected_sq_dists[k], abs_tol=1e-8)


@pytest.mark.parametrize(
    "plane_description, cell_params",
    [
        ("2D rectangular in ab-plane (c >> a, b)", (1.0, 1.4, 10.0, 90.0, 90.0, 90.0)),
        ("2D oblique in ab-plane (gamma = 72 deg)", (1.0, 1.3, 10.0, 90.0, 90.0, 72.0)),
        ("2D oblique in bc-plane (a >> b, c, alpha = 110 deg)", (10.0, 1.1, 1.5, 110.0, 90.0, 90.0)),
        ("2D oblique in ac-plane (b >> a, c, beta = 65 deg)", (1.2, 10.0, 0.9, 90.0, 65.0, 90.0)),
    ],
)
def test_2d_anisotropic_and_oblique_lattices(
    plane_description: str,
    cell_params: tuple[float, float, float, float, float, float],
):
    """Verify 2D lattices across different planes, side-length ratios, and oblique angles.

    System description:
        Quasi-2D crystals where one unit-cell axis is ``10.0`` Angstrom and the
        other two axes have unequal lengths and either orthogonal or oblique
        inter-axial angles (``alpha``, ``beta``, or ``gamma != 90``). Tested up
        to ``n = 8``.
    """
    unit_cell = UnitCell(*cell_params)
    site = Site(0.25, 0.4, 0.1, name=plane_description)
    n = 8

    shells = NeighbourSearch(site, [site], unit_cell, n=n,
                             max_iters=MAX_OUTER_LOOP_ITERATIONS).neighbour_sites()
    assert_shells_match_reference(shells, site, [site], unit_cell, n=n)


def test_2d_kagome_lattice_via_structure_and_spacegroup():
    """Verify nearest neighbours on a 2D Kagome lattice built with ``Structure`` and ``P6/mmm``.

    System description:
        A 2D Kagome net constructed using space group ``P6/mmm`` and
        ``UnitCell(1.0, 1.0, 10.0, gamma=120)`` from a single symmetry-generating
        site at Wyckoff position ``3f`` ``(0.5, 0.0, 0.0)``, which generates 3
        Kagome sites per unit cell: ``(0.5, 0, 0)``, ``(0, 0.5, 0)``, and
        ``(0.5, 0.5, 0)``.
        - When searching from one Kagome site across all 3 Kagome sites
          (``start_site`` in ``search_sites``), the shell sizes are:
          1st NN (``d = 0.5``): 4 sites; 2nd NN (``d = sqrt(3)/2``): 4 sites;
          3rd NN (``d = 1.0``): 6 sites.
        - When searching from the hexagon center ``(0, 0, 0)`` (NOT in
          ``search_sites``) to the 3 Kagome sites, the 1st shell is a 6-site
          hexagon at distance ``0.5``.
    """
    unit_cell = UnitCell(1.0, 1.0, 10.0, gamma=120.0)
    sg = get_spacegroup("P6/mmm")
    structure = Structure([Site(0.5, 0.0, 0.0, name="kagome")], unit_cell=unit_cell, spacegroup=sg)
    kagome_sites = structure.sites
    assert len(kagome_sites) == 3

    # 1. Start site is one of the Kagome sites (in search_sites)
    shells_in = NeighbourSearch(kagome_sites[0], kagome_sites, unit_cell, n=6,
                                max_iters=MAX_OUTER_LOOP_ITERATIONS).neighbour_sites()
    dists_in = assert_shells_match_reference(shells_in, kagome_sites[0], kagome_sites, unit_cell, n=6)
    assert [len(shells_in[k]) for k in range(4)] == [1, 4, 4, 6]
    assert math.isclose(dists_in[1], 0.5, abs_tol=1e-8)
    assert math.isclose(dists_in[2], math.sqrt(3) / 2.0, abs_tol=1e-8)
    assert math.isclose(dists_in[3], 1.0, abs_tol=1e-8)

    # 2. Start site is the hexagon center (0, 0, 0), NOT in search_sites
    hexagon_center = Site(0.0, 0.0, 0.0, name="center")
    shells_out = NeighbourSearch(hexagon_center, kagome_sites, unit_cell, n=5,
                                 max_iters=MAX_OUTER_LOOP_ITERATIONS).neighbour_sites()
    dists_out = assert_shells_match_reference(shells_out, hexagon_center, kagome_sites, unit_cell, n=5)
    assert len(shells_out[0]) == 0
    assert len(shells_out[1]) == 6
    assert math.isclose(dists_out[1], 0.5, abs_tol=1e-8)


# ===========================================================================
# 4. 3D Periodic Systems: Bravais Lattices, All Crystal Systems & Unit Cells
# ===========================================================================


@pytest.mark.parametrize(
    "lattice_name, sg_name, expected_sites_in_cell, expected_shell_sizes",
    [
        ("Simple Cubic (sc)", "Pm-3m", 1, [1, 6, 12, 8, 6, 24, 24]),
        ("Body-Centered Cubic (bcc)", "Im-3m", 2, [1, 8, 6, 12, 24, 8, 6]),
        ("Face-Centered Cubic (fcc)", "Fm-3m", 4, [1, 12, 6, 24, 12, 24, 8]),
    ],
)
def test_3d_cubic_bravais_lattices_via_spacegroup(
    lattice_name: str,
    sg_name: str,
    expected_sites_in_cell: int,
    expected_shell_sizes: list[int],
):
    """Verify analytical coordination shell sizes for 3D SC, BCC, and FCC crystals.

    System description:
        Cubic crystals (``a = b = c = 2.0`` Angstrom) generated from a single
        origin site ``(0, 0, 0)`` using space groups ``Pm-3m`` (1 site/cell),
        ``Im-3m`` (2 sites/cell), and ``Fm-3m`` (4 sites/cell). Checks the
        first ``n = 6`` coordination shells against known crystallographic
        coordination sequences.
    """
    unit_cell = UnitCell(2.0, 2.0, 2.0)
    sg = get_spacegroup(sg_name)
    structure = Structure([Site(0.0, 0.0, 0.0, name=lattice_name)], unit_cell=unit_cell, spacegroup=sg)
    sites = structure.sites
    assert len(sites) == expected_sites_in_cell

    n = 6
    shells = NeighbourSearch(sites[0], sites, unit_cell, n=n,
                             max_iters=MAX_OUTER_LOOP_ITERATIONS).neighbour_sites()
    assert_shells_match_reference(shells, sites[0], sites, unit_cell, n=n)
    assert [len(s) for s in shells] == expected_shell_sizes


@pytest.mark.parametrize(
    "system_label, cell_params",
    [
        ("Tetragonal (a=b!=c, 90/90/90)", (1.2, 1.2, 1.9, 90.0, 90.0, 90.0)),
        ("Orthorhombic (a!=b!=c, 90/90/90)", (1.0, 1.3, 1.7, 90.0, 90.0, 90.0)),
        ("Hexagonal (a=b!=c, 90/90/120)", (1.3, 1.3, 2.1, 90.0, 90.0, 120.0)),
        ("Rhombohedral acute (a=b=c, 60/60/60)", (1.5, 1.5, 1.5, 60.0, 60.0, 60.0)),
        ("Rhombohedral obtuse (a=b=c, 104/104/104)", (1.5, 1.5, 1.5, 104.0, 104.0, 104.0)),
        ("Monoclinic (a!=b!=c, beta=115)", (1.1, 1.4, 1.8, 90.0, 115.0, 90.0)),
        ("Triclinic acute (all lengths/angles distinct)", (1.1, 1.4, 1.8, 75.0, 82.0, 68.0)),
        ("Triclinic obtuse (all lengths/angles distinct)", (1.2, 1.5, 0.9, 105.0, 98.0, 112.0)),
        ("Triclinic mixed acute/obtuse", (1.0, 1.6, 1.25, 78.0, 108.0, 84.0)),
    ],
)
def test_3d_diverse_unit_cell_geometries(
    system_label: str,
    cell_params: tuple[float, float, float, float, float, float],
):
    """Verify ``NeighbourSearch`` across all 3D unit-cell shapes and inter-axial angles.

    System description:
        3D Bravais lattices spanning tetragonal, orthorhombic, hexagonal,
        acute and obtuse rhombohedral, monoclinic, and acute/obtuse/mixed
        triclinic unit cells with unequal side lengths ``(a, b, c)`` and
        non-90-degree angles ``(alpha, beta, gamma)``. Tested up to ``n = 8``
        both when ``start_site`` is in ``search_sites`` and when ``start_site``
        is a distinct interstitial site not in ``search_sites``.
    """
    unit_cell = UnitCell(*cell_params)
    lattice_site = Site(0.2, 0.35, 0.6, name=system_label)
    interstitial_site = Site(0.75, 0.15, 0.85, name="interstitial")
    n = 8

    # 1. start_site in search_sites
    shells_in = NeighbourSearch(lattice_site, [lattice_site], unit_cell, n=n,
                                max_iters=MAX_OUTER_LOOP_ITERATIONS).neighbour_sites()
    assert_shells_match_reference(shells_in, lattice_site, [lattice_site], unit_cell, n=n)

    # 2. start_site NOT in search_sites
    shells_out = NeighbourSearch(interstitial_site, [lattice_site], unit_cell, n=n,
                                 max_iters=MAX_OUTER_LOOP_ITERATIONS).neighbour_sites()
    assert_shells_match_reference(shells_out, interstitial_site, [lattice_site], unit_cell, n=n)


# ===========================================================================
# 5. High and Low Symmetry Multi-Site Structures (Small & Large Site Counts)
# ===========================================================================


def test_low_symmetry_p1_four_site_example():
    """Verify the 4-site low-symmetry (P1) example from ``nearest_neighbours.py``.

    System description:
        Reproduces the first ``__main__`` test case in ``nearest_neighbours.py``:
        ``UnitCell(1, 0.2, 1)`` with four low-symmetry sites at ``(0.2, 0.3, 0.4)``,
        ``(0.5, 0.5, 0.5)``, ``(0.2, 0.4, 0.6)``, and ``(0.8, 0.7, 0.6)``,
        searching from ``sites[1]`` up to ``n = 6``, plus searching from
        ``sites[1]`` into the other three sites (``start_site`` not in
        ``search_sites``).
    """
    unit_cell = UnitCell(1.0, 0.2, 1.0)
    sites = [
        Site(0.2, 0.3, 0.4, name="s0"),
        Site(0.5, 0.5, 0.5, name="s1"),
        Site(0.2, 0.4, 0.6, name="s2"),
        Site(0.8, 0.7, 0.6, name="s3"),
    ]
    n = 6

    # Start site in search_sites
    shells_in = NeighbourSearch(sites[1], sites, unit_cell, n=n,
                                max_iters=MAX_OUTER_LOOP_ITERATIONS).neighbour_sites()
    assert_shells_match_reference(shells_in, sites[1], sites, unit_cell, n=n)

    # Start site NOT in search_sites
    other_sites = [sites[0], sites[2], sites[3]]
    shells_out = NeighbourSearch(sites[1], other_sites, unit_cell, n=n,
                                 max_iters=MAX_OUTER_LOOP_ITERATIONS).neighbour_sites()
    assert_shells_match_reference(shells_out, sites[1], other_sites, unit_cell, n=n)


@pytest.mark.parametrize(
    "description, sg_name, cell_params, generator_coords, expected_site_count, n",
    [
        (
            "Low-symmetry triclinic P-1 with 2 general sites (4 sites/cell)",
            "P-1",
            (1.2, 1.5, 1.9, 78.0, 85.0, 108.0),
            [(0.18, 0.31, 0.42), (0.63, 0.12, 0.77)],
            4,
            8,
        ),
        (
            "Low-symmetry monoclinic P21/c general Wyckoff 4e (4 sites/cell)",
            "P21/c",
            (1.3, 1.7, 2.1, 90.0, 112.0, 90.0),
            [(0.19, 0.33, 0.27)],
            4,
            8,
        ),
        (
            "Medium-symmetry orthorhombic Pmmm general Wyckoff 8α (8 sites/cell)",
            "Pmmm",
            (1.4, 1.8, 2.3, 90.0, 90.0, 90.0),
            [(0.17, 0.23, 0.31)],
            8,
            7,
        ),
        (
            "High-symmetry rhombohedral R-3mR general Wyckoff 12i (12 sites/cell)",
            "R-3mR",
            (2.0, 2.0, 2.0, 70.0, 70.0, 70.0),
            [(0.14, 0.29, 0.47)],
            12,
            7,
        ),
        (
            "High-symmetry tetragonal P4/mmm general Wyckoff 16u (16 sites/cell)",
            "P4/mmm",
            (2.0, 2.0, 2.8, 90.0, 90.0, 90.0),
            [(0.13, 0.29, 0.37)],
            16,
            6,
        ),
        (
            "High-symmetry hexagonal P6/mmm general Wyckoff 24r (24 sites/cell)",
            "P6/mmm",
            (2.2, 2.2, 3.1, 90.0, 90.0, 120.0),
            [(0.15, 0.37, 0.23)],
            24,
            6,
        ),
        (
            "High-symmetry cubic Pm-3m general Wyckoff 48n (48 sites/cell)",
            "Pm-3m",
            (3.0, 3.0, 3.0, 90.0, 90.0, 90.0),
            [(0.12, 0.23, 0.37)],
            48,
            6,
        ),
    ],
)
def test_symmetric_structures_across_site_multiplicities(
    description: str,
    sg_name: str,
    cell_params: tuple[float, float, float, float, float, float],
    generator_coords: list[tuple[float, float, float]],
    expected_site_count: int,
    n: int,
):
    """Verify ``NeighbourSearch`` on ``Structure`` + ``SpaceGroup`` crystals from 4 to 48 sites/cell.

    System description:
        Crystals generated by applying space groups ranging from low symmetry
        (``P-1``, ``P21/c``) to high symmetry (``Pmmm``, ``R-3mR``, ``P4/mmm``,
        ``P6/mmm``, ``Pm-3m``) to general Wyckoff positions, producing between
        4 and 48 symmetry-equivalent sites per unit cell.
        Tests both:
        1. ``start_site`` in ``search_sites`` (searching all sites in the cell).
        2. ``start_site`` NOT in ``search_sites`` (searching ``sites[1:]`` from
           ``sites[0]``).
    """
    unit_cell = UnitCell(*cell_params)
    sg = get_spacegroup(sg_name)
    input_sites = [
        Site(x, y, z, name=f"gen_{idx}")
        for idx, (x, y, z) in enumerate(generator_coords)
    ]
    structure = Structure(input_sites, unit_cell=unit_cell, spacegroup=sg)
    all_sites = structure.sites
    assert len(all_sites) == expected_site_count, (
        f"Expected {expected_site_count} sites for {description}, got {len(all_sites)}"
    )

    start_site = all_sites[0]

    # 1. Start site in search_sites
    shells_in = NeighbourSearch(start_site, all_sites, unit_cell, n=n,
                                max_iters=MAX_OUTER_LOOP_ITERATIONS).neighbour_sites()
    assert_shells_match_reference(shells_in, start_site, all_sites, unit_cell, n=n)

    # 2. Start site NOT in search_sites
    other_sites = all_sites[1:]
    shells_out = NeighbourSearch(start_site, other_sites, unit_cell, n=n,
                                 max_iters=MAX_OUTER_LOOP_ITERATIONS).neighbour_sites()
    assert_shells_match_reference(shells_out, start_site, other_sites, unit_cell, n=n)


def test_start_site_near_cell_corner_finds_closer_neighbour_in_adjacent_cell():
    """Verify that a site near a unit-cell boundary finds its true nearest neighbour across the cell face/corner.

    System description:
        A triclinic unit cell ``UnitCell(1.5, 1.8, 2.1, 80, 95, 105)`` where
        ``start_site`` is near the ``(1, 1, 1)`` corner at fractional coordinates
        ``(0.92, 0.91, 0.93)`` and ``search_sites`` contains a single site near
        the origin at ``(0.05, 0.06, 0.04)``.
        The true 1st nearest neighbour of ``start_site`` is the image of
        ``search_sites[0]`` in unit cell ``(+1, +1, +1)`` (fractional separation
        ``(0.13, 0.15, 0.11)``), which is much closer than the image in cell
        ``(0, 0, 0)`` (fractional separation ``(-0.87, -0.85, -0.89)``).
    """
    unit_cell = UnitCell(1.5, 1.8, 2.1, alpha=80.0, beta=95.0, gamma=105.0)
    start_site = Site(0.92, 0.91, 0.93, name="corner_start")
    target_site = Site(0.05, 0.06, 0.04, name="origin_target")
    n = 6

    shells = NeighbourSearch(start_site, [target_site], unit_cell, n=n,
                             max_iters=MAX_OUTER_LOOP_ITERATIONS).neighbour_sites()
    assert_shells_match_reference(shells, start_site, [target_site], unit_cell, n=n)

    # The 1st nearest neighbour must be in the adjacent cell (+1, +1, +1)
    assert shells[1] == [(target_site, (1, 1, 1))]
