""" Sample deserialisation tests """
import numpy as np
import pytest

from pyspinw import *
from pyspinw.basis import angle_axis_rotation_matrix
from pyspinw.sample import Sample

site = LatticeSite(1/2,1/2,1/2, name="test_site")
structure = Structure([site], unit_cell=UnitCell(1,1,1))
hamiltonian = Hamiltonian(structure, [])

def test_single_crystal_serialisation():
    """ Single Crystal serialisation test """

    sample = SingleCrystal(hamiltonian)

    json = sample.serialise()

    deserialised = Sample.deserialise(json)

    assert isinstance(deserialised, SingleCrystal)

    # Basic check that the hamiltonian is what it should be
    assert hamiltonian.structure.sites[0].name == "test_site"


def test_twin_serialisation():
    """ Twin serialisation test """

    sample = Twin(hamiltonian, [1,1,1], 0.3)

    json = sample.serialise()

    deserialised = Sample.deserialise(json)

    assert isinstance(deserialised, Twin)

    # Multidomain checks
    assert np.all(np.array(sample.weights == np.array(deserialised.weights)))
    assert len(sample._transformations) == len(deserialised._transformations)
    for t1, t2 in zip(sample._transformations, deserialised._transformations):
        assert np.all(t1 == t2)

    ## check domains
    assert len(sample._domains) == len(deserialised._domains)
    for d1, d2 in zip(sample._domains, deserialised._domains):
        assert d1.weighting == d2.weighting
        assert np.all(d1.transformation == d2.transformation)

    # Twin specific checks
    assert np.all(sample._twinning_plane_normal == deserialised._twinning_plane_normal)
    assert sample._second_twin_fraction == deserialised._second_twin_fraction

    # Basic check that the hamiltonian is what it should be
    assert hamiltonian.structure.sites[0].name == "test_site"

def test_multidomain_serialisation():
    """ Multidomain serialisation test """

    domains = [
        CrystalDomain([[0,1,0],[1,0,0],[0,0,1]], 0.2),
        CrystalDomain(angle_axis_rotation_matrix(0.4, np.array([1,1,1], dtype=float)), 7),
        CrystalDomain(np.eye(3), np.sqrt(2))
    ]

    sample = Multidomain(hamiltonian, domains)

    json = sample.serialise()

    deserialised = Sample.deserialise(json)

    assert isinstance(deserialised, Multidomain)

    # Multidomain checks
    assert np.all(np.array(sample.weights == np.array(deserialised.weights)))
    assert len(sample._transformations) == len(deserialised._transformations)
    for t1, t2 in zip(sample._transformations, deserialised._transformations):
        assert np.all(t1 == t2)

    ## check domains
    assert len(sample._domains) == len(deserialised._domains)
    for d1, d2 in zip(sample._domains, deserialised._domains):
        assert d1.weighting == d2.weighting
        assert np.all(d1.transformation == d2.transformation)

    # Basic check that the hamiltonian is what it should be
    assert hamiltonian.structure.sites[0].name == "test_site"

def test_powder_serialisation():
    """ Powder serialisation test"""

    sample = Powder(hamiltonian)

    json = sample.serialise()

    deserialised = Sample.deserialise(json)

    assert isinstance(deserialised, Powder)

    # Basic check that the hamiltonian is what it should be
    assert hamiltonian.structure.sites[0].name == "test_site"