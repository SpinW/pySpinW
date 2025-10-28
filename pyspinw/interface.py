""" Helper functions for python interface """

from pyspinw.symmetry.group import database, NoSuchGroup, ExactMatch, PartialMatch


def spacegroup(search_string: str):

    if 'x' in search_string and 'y' in search_string and 'z' in search_string:
        # Not a spacegroup name, but might be a list of symmetry operations
        m = database.spacegroups_with_operations(search_string)

        if isinstance(m, ExactMatch):
            return m.spacegroup

        elif isinstance(m, PartialMatch):
            groups = m.spacegroups

            if len(groups) == 0:
                raise NoSuchGroup("Tried to find group by operations, no matching group found")

            elif len(groups) == 1:
                return groups[0]

            elif 1 < len(groups) <= 5:
                suggestions = ", ".join([group.symbol for group in groups[:-1]]) + " and " + groups[-1].symbol
                raise NoSuchGroup(f"Found multiple groups with these operations: {suggestions}")

            else:
                raise NoSuchGroup(f"Found {len(groups)} groups matching those operations.")

        else:
            raise ValueError("Expected `spacegroups_by_operations` to return ExactMatch or PartialMatch")

    else:
        # Try to get by name
        return database.spacegroup_by_name(search_string)


