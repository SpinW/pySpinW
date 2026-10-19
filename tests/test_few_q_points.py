""" The Python spin-wave path works with fewer q-points than worker processes (issue #253) """

import numpy as np
import pytest

import pyspinw.calculations.spinwave
from pyspinw.exchange import HeisenbergExchange
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.site import LatticeSite
from pyspinw.structure import Structure
from pyspinw.symmetry.unitcell import UnitCell


@pytest.mark.parametrize("n_q", [1, 2, 4, 5])
def test_python_path_with_fewer_q_points_than_processes(monkeypatch, n_q):
    """ With 20 cores the q-points are split into 5 chunks; fewer q-points than that used to fail """
    monkeypatch.setattr(pyspinw.calculations.spinwave.multiprocessing, "cpu_count", lambda: 20)

    site = LatticeSite(0, 0, 0, 0, 0, 1, name="S")
    hamiltonian = Hamiltonian(Structure([site], UnitCell(3.0, 8.0, 8.0)),
                              [HeisenbergExchange(site, site, j=-1.0, cell_offset=(1, 0, 0))])
    q_vectors = np.column_stack([np.linspace(0.1, 0.4, n_q), np.zeros(n_q), np.zeros(n_q)])

    energies, _ = hamiltonian._energies_and_intensities(q_vectors, use_rust=False)
    reference, _ = hamiltonian._energies_and_intensities(q_vectors, use_rust=True)

    assert np.array(energies).shape == (n_q, 2)
    assert np.sort(np.real(energies), axis=1) == pytest.approx(np.sort(np.real(reference), axis=1))
