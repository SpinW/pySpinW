import numpy as np
import pytest

from pyspinw import *

test_slices = [
    (1, 2, 3),
    (slice(None, None, None), 0, 0),
    (slice(None, None, None), 0, slice(2, None, None)),
    (slice(None, None, None), slice(None, None, None), slice(None, None, None)),
]

@pytest.mark.parametrize("input_type", [0,1,2,3], ids=["[full matrix]", "[tuple]", "[1D]", "[1x3]"])
@pytest.mark.parametrize("create_structure_first", [True, False], ids=["[structure first]", "[no structure]"])
@pytest.mark.parametrize("slices", test_slices, ids=["[tuple]", "[one slice]", "[two slices]", "[all slices]"])
def test_direct_supercell_manual_spin_setting(create_structure_first, input_type, slices):
    """ Checks the methods for setting spins """

    size = 2, 3, 5
    supercell = DirectSupercell(*size)

    def bool_indices(item):
        output = np.zeros(size, dtype=bool)
        output[item] = True
        return output.reshape(-1)


    site = LatticeSite(1/2,1/2,1/2,1,2,3)

    should_be_changed = bool_indices(slices)

    should_not_be_changed = ~should_be_changed

    n_changes = np.sum(should_be_changed)

    # Different ways of setting the input values
    match input_type:
        case 0:
            test_value = np.repeat(np.array([[4,5,6]]), n_changes, axis=0)
        case 1:
            test_value = 4,5,6
        case 2:
            test_value = np.array([4,5,6])
        case 3:
            test_value = np.array([[4,5,6]])
        case _:
            raise ValueError("Test values are 0,1,2,3")

    if create_structure_first:
        # Check that creating the structure (which coerces spin_data) works with it
        structure = Structure([site], UnitCell(1,1,1), supercell=supercell)
        supercell.spins_for(site).__setitem__(slices, test_value)

    else:
        supercell.spins_for(site).__setitem__(slices, test_value)

    assert np.all(site.spin_data[should_not_be_changed, 0] == 1)
    assert np.all(site.spin_data[should_not_be_changed, 1] == 2)
    assert np.all(site.spin_data[should_not_be_changed, 2] == 3)

    assert np.all(site.spin_data[should_be_changed, 0] == 4)
    assert np.all(site.spin_data[should_be_changed, 1] == 5)
    assert np.all(site.spin_data[should_be_changed, 2] == 6)


def test_direct_supercell_spin_setting_correctness_by_expanding():
    """ Checks that the calculations of the spin are correct """
    size = 2, 3, 5
    supercell = DirectSupercell(*size)

    site = LatticeSite(0,0,0)
    structure = Structure([site], UnitCell(1,1,1), supercell=supercell)

    for cell in supercell.cells():
        supercell.spins_for(site)[cell.as_tuple] = cell.as_tuple

    expanded = structure.expand()

    for site in expanded.sites:
        for i in range(3):
            position = site.ijk
            spin = site.spin
            assert np.isclose(size[i]*position[i], spin[i]), "Position in expanded cell should be spin / size"


rng = np.random.default_rng(2001)
def test_direct_supercell_spins_are_optimisable():
    """ Checks the calculations of derivatives for optimisation are correct """

    size = 2, 3, 5
    supercell = DirectSupercell(*size)

    site = LatticeSite(1/2,1/2,1/2)
    structure = Structure([site], UnitCell(1, 1, 1), supercell=supercell)

    supercell.spins_for(site)[:,:,:] = rng.normal(size=size + (3,)).reshape(-1, 3)

    hamiltonian = Hamiltonian(structure, [])

    field = [1,1,1]

    optimised = hamiltonian.ground_state(field=field)

    for site in optimised.expanded().structure.sites:
        spin_direction = site.spin_data[0,:] / np.sqrt(np.sum(site.spin_data[0,:]**2))

        # Unit vector dot [1,1,1] should be sqrt(3) for parallel, or -sqrt(3) for antiparallel
        assert np.isclose(np.dot(spin_direction, field), -np.sqrt(3)), \
            "Optimised spin should point in the opposite direction of the field"