# Docstring Overview Improvements Needed: `structure.py`

## `Structure._without_nonmagnetic_with_removed_uids` (line 365)
- **Current overview**: `Get a copy of the structure without magnetic sites, along with a list of uids for the removed sites.`
- **Issue**: The summary states "without magnetic sites", which is the opposite of the method's name and behavior (it removes *non-magnetic* sites whose spin is zero and keeps magnetic sites). It also says "list of uids" when the second return value is a `set[int]`.
- **Suggested clarification needed**: Change "without magnetic sites" to "without non-magnetic sites" and "list of uids" to "set of unique IDs".

## `Structure.sites_by_name` (line 417)
- **Current overview**: `Get sites where name matches regex.`
- **Issue**: Despite claiming to match a regular expression, the implementation performs a case-insensitive prefix check (`site.name.lower().startswith(name.lower())`).
- **Suggested clarification needed**: Clarify whether the docstring should describe a case-insensitive prefix match or whether the implementation was intended to use regex matching.

## `Structure.one_neighbour` (line 604)
- **Current overview**: `Get one, arbitrary, n-th nearest neighbours, along with cell offsets.`
- **Issue**: Plural "neighbours" and "cell offsets" in the summary, whereas the method returns a single `(Site, CellOffset)` tuple.
- **Suggested clarification needed**: Update the summary grammar to singular ("Get one arbitrary n-th nearest neighbour, along with its cell offset.").

## `Structure.updated` (line 699)
- **Current overview**: `Make a copy of this structure, but with selected (not None) fields replaced.`
- **Issue**: Accepts a `sites` parameter, but passes `self._input_sites` unconditionally to `Structure(...)`, so replacing `sites` has no effect on the returned structure's `_input_sites`.
- **Suggested clarification needed**: Clarify or fix how `sites` is handled in `Structure.updated`.
