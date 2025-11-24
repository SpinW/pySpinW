from pyspinw.site import LatticeSite
from pyspinw.symmetry.group import SpaceGroup, MagneticSpaceGroup
from pyspinw.symmetry.supercell import Supercell
from pyspinw.symmetry.unitcell import UnitCell



class Structure:
    def __init__(self,
                 sites: list[LatticeSite],
                 spacegroup: SpaceGroup | MagneticSpaceGroup,
                 unit_cell: UnitCell,
                 supercell: Supercell):

        self._input_sites = sites
        self._spacegroup = spacegroup
        self._unit_cell = unit_cell
        self._supercell = supercell

        self._sites: list[LatticeSite] = self._build_sites()

    def _extended_sites(self)  -> list[LatticeSite]:
        pass

    def _build_sites(self):
        self._sites = self._extended_sites()

    @property
    def sites(self) -> list[LatticeSite]:
        return self._sites.copy()

    @sites.setter
    def sites(self, sites):
        self._input_sites = sites
        self._sites = self._build_sites()

    @property
    def spacegroup(self) -> SpaceGroup | MagneticSpaceGroup:
        return self._spacegroup

    @spacegroup.setter
    def spacegroup(self, spacegroup: SpaceGroup | MagneticSpaceGroup):
        self._spacegroup = spacegroup
        self._build_sites()

    @property
    def unit_cell(self) -> UnitCell:
        return self._unit_cell

    @unit_cell.setter
    def unit_cell(self, unit_cell: UnitCell):
        self._unit_cell = unit_cell
        self._build_sites()

    @property
    def supercell(self) -> Supercell:
        return self._supercell

    @supercell.setter
    def supercell(self, supercell: Supercell):
        self._supercell = supercell
        self._build_sites()