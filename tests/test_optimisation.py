import pytest

import numpy as np


from pyspinw.calculations.optimisation.energy_minimisation import ClassicalEnergyMinimisation, Free, Fixed, Planar
from pyspinw.interface import generate_exchanges, axis_anisotropies
from pyspinw.coupling import HeisenbergCoupling
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.site import LatticeSite
from pyspinw.structures import Structure
from pyspinw.symmetry.supercell import TrivialSupercell, SummationSupercell, PropagationVector, \
    CommensuratePropagationVector
from pyspinw.symmetry.unitcell import UnitCell

@pytest.mark.parametrize("size", [0.1, 0.4])
def test_jitter_free(size):
    """ Check that the jitter method does jittering as expected for free sites"""
    rng = np.random.default_rng(911)
    start_moments = rng.normal(0,1,(30, 1, 3))
    start_moments /= np.sqrt(np.sum(start_moments**2, axis=2)).reshape(-1, 1, 1)

    assert start_moments.shape == (30, 1, 3), "Input shape must be n_sites-by-n_components-by-3"

    sites = [LatticeSite(m1, m2, m3, m1, m2, m3) for m1, m2, m3 in start_moments.reshape(30, 3)]

    unit_cell = UnitCell(1, 1, 1, gamma=60)
    s = Structure(sites, unit_cell=unit_cell, supercell=TrivialSupercell())

    hamiltonian = Hamiltonian(s, [])

    hamiltonian.print_summary()

    minimiser = ClassicalEnergyMinimisation(hamiltonian, constraints=Free, seed=91210)

    minimiser.jitter(size)

    amounts = np.sum(minimiser.moment_data * start_moments, axis=2)

    assert np.allclose(amounts, np.cos(size))

@pytest.mark.parametrize("size", [0.1, 0.4])
def test_jitter_planar(size):
    rng = np.random.default_rng(911)
    start_moments = rng.normal(0,1,(30, 1, 3))
    start_moments /= np.sqrt(np.sum(start_moments**2, axis=2)).reshape(-1, 1, 1)

    assert start_moments.shape == (30, 1, 3), "Input shape must be n_sites-by-n_components-by-3"

    sites = [LatticeSite(m1, m2, m3, m1, m2, m3) for m1, m2, m3 in start_moments.reshape(30, 3)]

    unit_cell = UnitCell(1, 1, 1, gamma=60)
    s = Structure(sites, unit_cell=unit_cell, supercell=TrivialSupercell())

    hamiltonian = Hamiltonian(s, [])

    hamiltonian.print_summary()

    minimiser = ClassicalEnergyMinimisation(hamiltonian, constraints=Planar([0,0,1]), seed=91210)

    minimiser.jitter(size)

    # Project to check angles
    start_project = start_moments[:, :, :2]
    start_project /= np.sqrt(np.sum(start_project**2, axis=2)).reshape(-1, 1, 1)

    jittered_project = minimiser.moment_data[:, :, :2]
    jittered_project /= np.sqrt(np.sum(jittered_project**2, axis=2)).reshape(-1, 1, 1)

    amounts = np.sum(start_project * jittered_project, axis=2)

    assert np.allclose(amounts, np.cos(size))



@pytest.mark.parametrize("dim", [2, 3])
@pytest.mark.parametrize("shape", [(1,), (5,), (3, 4), (1, 2), (6, 1), (2, 2, 3)])
def test_random_orientations(shape, dim):

    # We want to check that the safety thing is only replacing rows that are all zeros
    # the plan then is return an integer that is different each time it is called, this
    # way we track which rows have been replaced

    class RNGMock:
        def __init__(self):
            self.counter = 1.0
            self.rng = np.random.default_rng(1979)

        def normal(self, mu, sigma, shape):
            out = np.array(self.rng.integers(0, 2, size=shape), dtype=float)
            out *= self.counter
            self.counter += 1
            return out

    rng = RNGMock()

    random_vectors = ClassicalEnergyMinimisation._random_orientations(rng, shape, dim)

    assert random_vectors.shape == shape + (dim, ), "Output array should have requested dimensions"

    def check_row(row):
        zeros = row == 0
        assert sum(zeros) != len(row), "Row should not be all zeros"

        others = row[~zeros]

        # others[0] should exist if the previous assertion has not failed
        assert np.all(others[0] == others[1:]), "All components should be the same"

        assert np.isclose(np.sum(row**2), 1), "Length should be 1"

    if len(shape) == 1:
        for i in range(shape[0]):
            check_row(random_vectors[i, :])

    elif len(shape) == 2:
        for i in range(shape[0]):
            for j in range(shape[1]):
                check_row(random_vectors[i,j,:])

    elif len(shape) == 3:
        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    check_row(random_vectors[i,j,k,:])

    else:
        assert False, "test setup is bad"

def test_randomisation_free():
    """ Check that the randomisation method randomises as expected for free sites"""
    n_sites = 30

    rng = np.random.default_rng(911)
    start_moments = rng.normal(0, 1, (n_sites, 2, 3))

    assert start_moments.shape == (30, 2, 3), "Input shape must be n_sites-by-n_components-by-3"

    sites = [LatticeSite(x[0,0], x[0,1], x[0,2], supercell_moments=x) for x in [
                start_moments[i, :, :] for i in range(n_sites) ]]

    unit_cell = UnitCell(1, 1, 1, gamma=60)
    s = Structure(sites, unit_cell=unit_cell, supercell=SummationSupercell([
            CommensuratePropagationVector(1,1,1),
            CommensuratePropagationVector(1,1,1)]))

    hamiltonian = Hamiltonian(s, [])

    hamiltonian.print_summary()

    minimiser = ClassicalEnergyMinimisation(hamiltonian, constraints=Free, seed=91210)

    minimiser.randomise()

    start_moment_square_mags = np.sum(start_moments**2, axis=2)

    randomised_moment_square_mags = np.sum(minimiser.moment_data**2, axis=2)

    assert np.allclose(start_moment_square_mags, randomised_moment_square_mags)

@pytest.mark.parametrize("axis", [[1,1,1], [1,2,3], [0,1,0], [0,0,1]])
def test_randomisation_planar(axis):
    """ Check that the randomisation method randomises as expected for free sites"""
    n_sites = 30

    rng = np.random.default_rng(911)
    start_moments = rng.normal(0, 1, (n_sites, 2, 3))

    assert start_moments.shape == (30, 2, 3), "Input shape must be n_sites-by-n_components-by-3"

    sites = [LatticeSite(x[0,0], x[0,1], x[0,2], supercell_moments=x) for x in [
                start_moments[i, :, :] for i in range(n_sites) ]]

    unit_cell = UnitCell(1, 1, 1, gamma=60)
    s = Structure(sites, unit_cell=unit_cell, supercell=SummationSupercell([
            CommensuratePropagationVector(1,1,1),
            CommensuratePropagationVector(1,1,1)]))

    hamiltonian = Hamiltonian(s, [])

    hamiltonian.print_summary()

    constraint = Planar(axis)
    minimiser = ClassicalEnergyMinimisation(hamiltonian, constraints=constraint, seed=91210)

    minimiser.randomise()

    start_moment_square_mags = np.sum(start_moments**2, axis=2)

    randomised_moment_square_mags = np.sum(minimiser.moment_data**2, axis=2)

    assert np.allclose(start_moment_square_mags, randomised_moment_square_mags), "Magnitudes should be conserved"

    prod = minimiser.moment_data * np.array(axis).reshape((1, 1, 3))
    moment_dot_axis = np.sum(prod, axis=2)

    assert np.allclose(moment_dot_axis, 0), "Moment should be perpendicular to specified axis"


@pytest.mark.parametrize("moment_size", [0.5, 1.0, 1.5])
def test_simple_ferromagnet_fixed_free(moment_size):
    unit_cell = UnitCell(1, 1, 1, gamma=60)

    x1 = LatticeSite(0, 0, 0.5, 0, 0, moment_size, name="X1")
    x2 = LatticeSite(0, 0, 0, 0, moment_size, 0, name="X2")

    sites = [x1, x2]

    s = Structure(sites, unit_cell=unit_cell, supercell=TrivialSupercell())

    exchanges = generate_exchanges(sites=sites,
                          unit_cell=unit_cell,
                          max_distance=0.6,
                          coupling_type=HeisenbergCoupling,
                          j=-1)

    hamiltonian = Hamiltonian(s, exchanges)

    hamiltonian.print_summary()

    minimiser = ClassicalEnergyMinimisation(hamiltonian, constraints=[Fixed, Free], field=np.array([0.0, 0.0, 0.0]))

    minimiser.minimise(verbose=True)

    assert np.allclose(minimiser.moment_data[1, 0, :], [0, 0, moment_size], atol=1e-4)



def test_simple_antiferromagnet_fixed_planar():
    unit_cell = UnitCell(1, 1, 1, gamma=60)

    x1 = LatticeSite(0, 0, 0.5, 0, 0, 1, name="X1")
    x2 = LatticeSite(0, 0, 0, 0, 1, 0, name="X2")

    sites = [x1, x2]

    s = Structure(sites, unit_cell=unit_cell, supercell=TrivialSupercell())

    exchanges = generate_exchanges(sites=sites,
                          unit_cell=unit_cell,
                          max_distance=0.6,
                          coupling_type=HeisenbergCoupling,
                          j=1)

    hamiltonian = Hamiltonian(s, exchanges)

    hamiltonian.print_summary()

    minimiser = ClassicalEnergyMinimisation(hamiltonian, constraints=[Fixed, Planar([1,0,0])], field=np.array([0.0, 0.0, 0.0]))

    minimiser.minimise(verbose=True)

    assert np.allclose(minimiser.moment_data[1, :], [0, 0, -1], atol=1e-5)

@pytest.mark.parametrize("moment_size", [0.5, 1.0, 1.5])
@pytest.mark.parametrize("axis", [[1,0,0],[1,1,1]])
def test_simple_antiferromagnet_free_planar(axis, moment_size):

    axis = np.array(axis, dtype=float)
    axis /= np.sqrt(np.sum(axis**2))

    unit_cell = UnitCell(1, 1, 1, gamma=60)

    x1 = LatticeSite(0, 0, 0.5, 0, 0, moment_size, name="X1")
    x2 = LatticeSite(0, 0, 0, 0, moment_size, 0, name="X2")

    sites = [x1, x2]

    s = Structure(sites, unit_cell=unit_cell, supercell=TrivialSupercell())

    exchanges = generate_exchanges(sites=sites,
                          unit_cell=unit_cell,
                          max_distance=0.6,
                          coupling_type=HeisenbergCoupling,
                          j=1)

    hamiltonian = Hamiltonian(s, exchanges)

    hamiltonian.print_summary()

    minimiser = ClassicalEnergyMinimisation(hamiltonian, constraints=[Free, Planar(axis)], field=np.array([0.0, 0.0, 0.0]))

    minimiser.minimise(verbose=True)


    assert np.isclose(np.dot(minimiser.moment_data[0, :], minimiser.moment_data[1, :]), -1)

    # We would like to check something like this - but its not true if they don't start off perpendicular to axis:
    #   assert np.isclose(np.dot(axis, minimiser.moments[1, :]), 0)
    # Instead check component is the same as when it started

    assert np.isclose(np.dot(axis, minimiser.moment_data[1, :]), np.dot(axis, sites[1].base_moment))


@pytest.mark.parametrize("axis", [[1,0,0],[1,1,1],[0,1,0],[-1,-1,0]])
@pytest.mark.parametrize("a", [-1, 1])
def test_anisotropies_free(axis, a):

    axis = np.array(axis, dtype=float)
    axis /= np.sqrt(np.sum(axis**2))

    unit_cell = UnitCell(1, 1, 1, gamma=60)

    x = LatticeSite(0, 0, 0.5, 0, 0, 1, name="X")

    sites = [x]

    s = Structure(sites, unit_cell=unit_cell, supercell=TrivialSupercell())

    ai = axis_anisotropies(sites, a, axis)

    hamiltonian = Hamiltonian(s, couplings=[], anisotropies=ai)

    hamiltonian.print_summary()

    minimiser = ClassicalEnergyMinimisation(hamiltonian, constraints=[Free], field=np.array([0.0, 0.0, 0.0]))

    minimiser.minimise(verbose=True)

    if a < 0:
        assert np.isclose(np.abs(np.dot(axis, minimiser.moment_data[0, :])), 1) # Should be aligned with axis (+-1)
    elif a > 0:
        assert np.isclose(np.dot(axis, minimiser.moment_data[0, :]), 0, atol=1e-4) # Should be in plane perpendicular to axis
    else:
        pass


@pytest.mark.parametrize("axis", [[1,0,0],[1,1,1],[0,1,0],[-1,-1,0]])
def test_anisotropies_planar(axis):
    """ Fixed to plane which it starts in

    if anisotropy forces a plane, it should always find it,
    if not, it won't be aligned
    """
    axis = np.array(axis, dtype=float)
    axis /= np.sqrt(np.sum(axis**2))

    unit_cell = UnitCell(1, 1, 1, gamma=60)

    x = LatticeSite(0, 0, 0.5, 0, 0, 1, name="X")

    sites = [x]

    s = Structure(sites, unit_cell=unit_cell, supercell=TrivialSupercell())

    ai = axis_anisotropies(sites, 1, axis)

    hamiltonian = Hamiltonian(s, couplings=[], anisotropies=ai)

    hamiltonian.print_summary()

    minimiser = ClassicalEnergyMinimisation(hamiltonian, constraints=[Planar([1,0,0])], field=np.array([0.0, 0.0, 0.0]))

    minimiser.minimise(verbose=True)

    assert np.isclose(np.dot(axis, minimiser.moment_data[0, :]), 0, atol=1e-4) # Should be in plane perpendicular to axis

rng = np.random.default_rng(118)
test_fields = rng.normal(0,1, (20, 3))

@pytest.mark.parametrize("field", test_fields)
def test_field_free(field):

    unit_cell = UnitCell(1, 1, 1, gamma=60)

    x1 = LatticeSite(0.1, 0, 0.5, 0, 1, 0, name="X1")
    x2 = LatticeSite(0.2, 0, 0.5, 1, 0, 0, name="X2")
    x3 = LatticeSite(0.3, 0, 0.5, 1, 1, 0, name="X3")
    x4 = LatticeSite(0.4, 0, 0.5, 1, -1, 0, name="X4")
    x5 = LatticeSite(0.5, 0, 0.5, -1, -1, 0, name="X5")

    sites = [x1, x2, x3, x4, x5]

    s = Structure(sites, unit_cell=unit_cell, supercell=TrivialSupercell())

    hamiltonian = Hamiltonian(s, couplings=[])

    hamiltonian.print_summary()

    minimiser = ClassicalEnergyMinimisation(hamiltonian, constraints=Free,
                                            field=field)

    minimiser.minimise(verbose=True)


    # Normalise moments and field
    field_norm = field / np.sqrt(np.sum(field**2))
    moments_norm = minimiser.moment_data.copy()
    moments_norm /= np.sqrt(np.sum(moments_norm**2, axis=1)).reshape(-1, 1)

    assert np.allclose(field_norm + moments_norm, 0.0, atol=1e-4)


@pytest.mark.parametrize("field", test_fields)
def test_field_planar(field):

    unit_cell = UnitCell(1, 1, 1, gamma=60)

    x1 = LatticeSite(0.1, 0, 0.5, 0, 1, 0, name="X1")
    x2 = LatticeSite(0.2, 0, 0.5, 1, 0, 0, name="X2")
    x3 = LatticeSite(0.3, 0, 0.5, 1, 1, 0, name="X3")
    x4 = LatticeSite(0.4, 0, 0.5, 1, -1, 0, name="X4")
    x5 = LatticeSite(0.5, 0, 0.5, -1, -1, 0, name="X5")

    sites = [x1, x2, x3, x4, x5]

    s = Structure(sites, unit_cell=unit_cell, supercell=TrivialSupercell())

    hamiltonian = Hamiltonian(s, couplings=[])

    hamiltonian.print_summary()

    minimiser = ClassicalEnergyMinimisation(hamiltonian, constraints=Planar([0, 0, 1]),
                                            field=field)

    minimiser.minimise(verbose=True)


    # Project to check angles
    field_project = field[:2]
    field_project /= np.sqrt(np.sum(field_project**2))
    field_project = field_project.reshape(1, 2)

    jittered_project = minimiser.moment_data[:, :2]
    jittered_project /= np.sqrt(np.sum(jittered_project**2, axis=1)).reshape(-1, 1)

    dot_product = np.sum(field_project * jittered_project, axis=1)

    assert np.allclose(dot_product, -1)