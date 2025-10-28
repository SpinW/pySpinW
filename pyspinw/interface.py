""" Helper functions for python interface """

from pyspinw.symmetry.group import database, NoSuchGroup, ExactMatch, PartialMatch


def spacegroup(search_string: str):
    """ Get a spacegroup by name or symmetry operations

    The searches are whitespace and case insensitive.

    Examples:
        The following are equivalent:
            spacegroup("P1")
            spacegroup("p1")
            spacegroup("x,y,z")

        The following are equivalent:
            spacegroup("P-1")
            spacegroup("p-1")
            spacegroup("x,y,z; -x,-y,-z")

        Spacegroups with multiple settings sometimes need to have it specified, but some don't:
            spacegroup("R3")  **fails**
            spacegroup("R3H")  **hexagonal setting**
            spacegroup("R3R")  **rhombohedral setting**

            spacegroup("B2/m") **this setting of C2/m**
            spacegroup("B 1 2/m 1") **another way for the same group**

            spacegroup("P 4/n 2/b 2/m : 1") **setting 1**
            spacegroup("P 4/n 2/b 2/m : 1") **setting 2**

    """
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


