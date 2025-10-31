""" Conversions and conventions for working with spacegroup string data"""

import re

import spglib

def canonise_string(name: str):
    """ Make name into a form for searching """
    # remove spaces
    name = re.sub(r"\s+", "", name)

    # make lowercase
    return name.lower()

def spacegroup_conventions():
    """ Get naming data from the spglib database that can be used for looking up spacegroups."""
    group_name_to_index = {}

    for i in range(1, 531):

        group_data = spglib.get_spacegroup_type(i)

        # Take the international name data and split by equality

        names = [s.strip() for s in group_data.international.split("=")]

        # Check that it is correct to use the first name in the lookup
        assert canonise_string(names[0]) == canonise_string(group_data.international_short)


        # Does the name contain the setting choice? this depends on what kind of choice it is
        # If the choice is 1/2, or R/H it doesn't

        # settings 1 and 2
        if group_data.choice.startswith("1") or group_data.choice.startswith("2"):
            names = [name + " : " + group_data.choice[:1] for name in names]

        # Rhombohedral/Hexagonal settings
        if group_data.choice in ["H", "R"]:
            names = [name + " " + group_data.choice for name in names]

        # Use the first entry as a short name if it is not there already
        if names[0] not in group_name_to_index:
            group_name_to_index[names[0]] = i


        # The rest we can add safely
        for name in names[1:]:
            group_name_to_index[name] = i

    for i in range(1, 531):
        group_data = spglib.get_spacegroup_type(i)

        # This is a little bit tricky, we need to make sure that the full name isn't already there
        # It might be there already, and point to a different group (because the international full
        # doesn't include the setting)
        full_name = group_data.international_full

        if group_data.choice.startswith("1") or group_data.choice.startswith("2"):
            full_name += " : " + group_data.choice[:1]

        if group_data.choice in ["R", "H"]:
            full_name += " " + group_data.choice

        if full_name not in group_name_to_index:
            group_name_to_index[full_name] = i

    # Make a lookup for the cannonical names, should be a one-to-one mapping
    canonical_name_to_index: dict[str, int] = {}
    canonical_name_to_group_name: dict[str, str] = {}
    for name, index in group_name_to_index.items():
        canonical_name = canonise_string(name)

        # check for collisions
        assert canonical_name not in canonical_name_to_index, "Conflicting canonical name"

        canonical_name_to_index[canonical_name] = index

        canonical_name_to_group_name[canonical_name] = name

    return canonical_name_to_group_name, canonical_name_to_index

if __name__ == "__main__":
    canonical_names_to_group_name, canonical_names_to_index = spacegroup_conventions()

    for key in canonical_names_to_group_name:
        print(key, "->", canonical_names_to_group_name[key])

