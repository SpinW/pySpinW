""" Tests for serialisation of the Hamiltonian class """
import numpy as np
import pytest

from pyspinw import *

def test_hamiltonian_serialisation():

    sites = [
        LatticeSite(0.1, 0.2, 0.3, 0, 0, 1, name="A"),
        LatticeSite(0.5, 0.5, 0, 1, 0, 0, name="B")
    ]

    sg = spacegroup("P4/mmm")
    unit_cell = UnitCell(1,1,1)
    supercell = TiledSupercell(scaling=(2,3,5))

    structure = Structure(sites, unit_cell=unit_cell, supercell=supercell, spacegroup=sg)

    anisotropies = [AxisMagnitudeAnisotropy(sites[0], direction=(1,0,0), a=-1),
                    Anisotropy(sites[1], anisotropy_matrix=np.eye(3))]

    exchanges = [HeisenbergExchange(sites[0], sites[0], cell_offset=(0,2,2), j=1),
                 HeisenbergExchange(sites[1], sites[0], j=1)]

    hamiltonian = Hamiltonian(structure, exchanges, anisotropies).symmetry_filled()

    json = hamiltonian.serialise()

    deserialised = Hamiltonian.deserialise(json)

    assert isinstance(deserialised, Hamiltonian)

    # We can assume the components are individually serialised correctly,
    # but we need to check exchanges and anisotropies are correctly assigned

    assert len(hamiltonian.structure.sites) == len(deserialised.structure.sites)

    # Exchange setup

    # need to check that things are populated, otherwise zip will have length zero and no tests will happen
    assert len(hamiltonian.exchanges) == len(deserialised.exchanges)

    for original, copied in zip(hamiltonian.exchanges, deserialised.exchanges):

        # Site 1
        assert np.all(original.site_1.ijk == copied.site_1.ijk)
        assert original.site_1.name == copied.site_1.name

        # Site 2
        assert np.all(original.site_2.ijk == copied.site_2.ijk)
        assert original.site_2.name == copied.site_2.name

    # Anisotropy setup

    # need to check that things are populated, otherwise zip will have length zero and no tests will happen
    assert len(hamiltonian.anisotropies) == len(deserialised.anisotropies)

    for original, copied in zip(hamiltonian.anisotropies, deserialised.anisotropies):
        assert np.all(original.site.ijk == copied.site.ijk)
        assert original.site.name == copied.site.name

        assert original.__class__ == copied.__class__


