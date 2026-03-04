""" Tests for expansion routines """
import numpy as np

from pyspinw.anisotropy import AxisMagnitudeAnisotropy
from pyspinw.coupling import HeisenbergCoupling
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.site import LatticeSite
from pyspinw.structures import Structure
from pyspinw.symmetry.group import database
from pyspinw.symmetry.supercell import TrivialSupercell
from pyspinw.symmetry.unitcell import UnitCell


def test_expansion_chain_no_rot():
    """ Test the expansion with a simple chain and simple duplication """

    site = LatticeSite(0,0,0,0,0,1)

    unexpanded_unit_cell = UnitCell(1,1,1)
    unexpanded_sites = [site]
    unexpanded_couplings = [
        HeisenbergCoupling(site, site, cell_offset=(1,0,0), j=-1, name="X"),
        HeisenbergCoupling(site, site, cell_offset=(0,1,0), j=-1, name="Y"),
        HeisenbergCoupling(site, site, cell_offset=(0,0,1), j=-1, name="Z")]

    unexpanded_anisotropies = [AxisMagnitudeAnisotropy(site, a=3)]

    unexpanded_structure = Structure(
                            sites=unexpanded_sites,
                            unit_cell=unexpanded_unit_cell,
                            spacegroup=database.spacegroup_by_name("P1"),
                            supercell=TrivialSupercell((1,2,3)))

    unexpanded_hamiltonian = Hamiltonian(
        structure=unexpanded_structure,
        couplings=unexpanded_couplings,
        anisotropies=unexpanded_anisotropies)

    expanded_hamiltonian = unexpanded_hamiltonian.expanded()

    # Check the unit cell

    assert expanded_hamiltonian.structure.unit_cell.a == 1
    assert expanded_hamiltonian.structure.unit_cell.b == 2
    assert expanded_hamiltonian.structure.unit_cell.c == 3

    # Check the sites
    ## First that there is the right number of them
    assert len(expanded_hamiltonian.structure.sites) == 1 * 2 * 3

    ## That they have the right positions
    for site in expanded_hamiltonian.structure.sites:
        ### Moment
        assert np.all(site.moment_data == np.array([0,0,1], dtype=float))

        ### Position, can't assume they're ordered, but they should be in one of a few positions
        assert site.i % 1 == 0        # Integer
        assert (site.j * 2) % 1 == 0  # Multiple of 1/2
        assert (site.k * 3) % 1 == 0  # Multiple of 1/3

        ### Check that the position is within [0, 1)
        assert 0 <= site.i < 1
        assert 0 <= site.j < 1
        assert 0 <= site.k < 1

    ## That they're all different
    for i, site_1 in enumerate(expanded_hamiltonian.structure.sites):
        for site_2 in expanded_hamiltonian.structure.sites[:i]:
            assert not (site_1.i == site_2.i and site_1.j == site_2.j and site_1.k == site_2.k)

    # Check the anisotropies have different sites
    for i, anisotropy_1 in enumerate(expanded_hamiltonian.anisotropies):
        for anisotropy_2 in expanded_hamiltonian.anisotropies[:i]:
            assert anisotropy_1.site._unique_id != anisotropy_2.site._unique_id

    # Check number of couplings
    assert len(expanded_hamiltonian.couplings) == len(unexpanded_hamiltonian.couplings) * 1 * 2 * 3

    # Check the couplings have the right intersite vector
    for coupling in expanded_hamiltonian.couplings:
        match coupling.name:
            case "X":
                expected_offset = np.array([1, 0, 0], dtype=float)
            case "Y":
                expected_offset = np.array([0, 1, 0], dtype=float)
            case "Z":
                expected_offset = np.array([0, 0, 1], dtype=float)
            case _:
                raise ValueError("Expected names 'X', 'Y' or 'Z'")

        position_1 = expanded_hamiltonian.structure.unit_cell.fractional_to_cartesian(coupling.site_1.ijk)
        position_2 = expanded_hamiltonian.structure.unit_cell.fractional_to_cartesian(coupling.site_2.ijk)

        offset = expanded_hamiltonian.structure.unit_cell.fractional_to_cartesian(coupling.cell_offset.vector)

        assert np.all(position_2 + offset - position_1 == expected_offset)



