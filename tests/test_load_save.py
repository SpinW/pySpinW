import pytest
from pyspinw import load, save, Site, Structure, UnitCell, Hamiltonian, SingleCrystal, spacegroup, load_spw

@pytest.mark.parametrize("add_extension", [True, False], ids=["with extension", "without extension"])
def test_load_save(tmp_path, add_extension):
    """ Check that loading and saving of internal structures works """

    site = Site(1 / 2, 1 / 2, 1 / 2, name="test_site")
    structure = Structure([site], unit_cell=UnitCell(1, 1, 2), spacegroup=spacegroup("p4mm"))
    hamiltonian = Hamiltonian(structure, [])
    sample = SingleCrystal(hamiltonian)

    filename_input = f"test_{add_extension}"

    save(sample, tmp_path / filename_input, add_extension=add_extension)
    if add_extension:
        to_check = load(tmp_path / f"{filename_input}.spw")
    else:
        to_check = load_spw(tmp_path / filename_input)


    # Do the check, mostly covered by specific serialisation things
    assert sample.__class__ == to_check.__class__
    assert structure.spacegroup.name == to_check.hamiltonian.structure.spacegroup.name
    assert site.name == to_check.hamiltonian.structure.sites[0].name


def test_load_save_alt_extension(tmp_path):
    """ Check that loading and saving of internal structures works """

    site = Site(1 / 2, 1 / 2, 1 / 2, name="test_site")
    structure = Structure([site], unit_cell=UnitCell(1, 1, 2), spacegroup=spacegroup("p4mm"))
    hamiltonian = Hamiltonian(structure, [])
    sample = SingleCrystal(hamiltonian)

    filename_input = f"test.json"

    save(sample, tmp_path / filename_input)
    to_check = load_spw(tmp_path / filename_input)

    # Do the check, mostly covered by specific serialisation things

    assert sample.__class__ == to_check.__class__
    assert structure.spacegroup.name == to_check.hamiltonian.structure.spacegroup.name
    assert site.name == to_check.hamiltonian.structure.sites[0].name




