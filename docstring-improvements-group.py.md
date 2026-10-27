# Docstring Overview Improvements Needed: `group.py`

## `SymmetryGroup.implied_sites_for` (line 38)
- **Current overview**: `Find "duplicate" sites of a given site.`
- **Issue**: Notice that `SymmetryGroup.implied_sites_for` accesses `self.operations`, which is not declared on `SymmetryGroup` (only `space_operations` is declared as an abstract property on `SymmetryGroup`), and "duplicate" sites actually means symmetry-equivalent sites within the unit cell.
- **Suggested clarification needed**: Clarify that this method generates all distinct symmetry-equivalent `ImpliedLatticeSite` instances within the unit cell for a given parent `site`.

## `SpacegroupDatabase.all_spacegroup_settings` (line 583)
- **Current overview**: `Get a list of all spacegroups that correspond to a given setting.`
- **Issue**: The parameter `number` is the International Tables space group number (1 to 230), and the method returns all settings (Hall numbers) corresponding to that space group number, not space groups corresponding to a setting.
- **Suggested clarification needed**: Update the overview to state that it returns all `SpaceGroup` setting variants corresponding to a given International Tables space group number.
