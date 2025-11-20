""" Loads data for lookups of spacegroup """

from importlib import resources

def load_formatted_aliases() -> dict[int, list[str]]:
    """ Load the aliases for all the spacegroup settings """

    output = {}

    with resources.open_text(
            "pyspinw.symmetry.data",
            "spacegroup_aliases_formatted.txt") as file:


        for line in file:
            parts = line.split(";")
            number = int(parts[0])
            aliases = [part.strip() for part in parts[2:]]

            output[number] = aliases


def load_canonical_aliases() -> dict[str, int]:
    """ Load the aliases for all the spacegroup settings in the form used for searching """

    output = {}

    with resources.open_text(
            "pyspinw.symmetry.data",
            "spacegroup_aliases_canonised.txt") as file:


        for line in file:
            parts = line.split(";")
            number = int(parts[0])
            aliases = [part.strip() for part in parts[2:]]

            for alias in aliases:
                output[alias] = number

    return output

def load_canonised_to_formatted() -> dict[str, str]:
    """ Load the mapping to convert lookup string to a nicely formatted string """
    output = {}

    with resources.open_text(
            "pyspinw.symmetry.data",
            "spacegroup_canonised_to_formatted.txt") as file:
        for line in file:
            parts = line.split(";")
            canonised = parts[0].strip()
            formatted = parts[1].strip()

            output[canonised] = formatted

    return output

def load_preferred_names() -> dict[int, str]:
    """ Load our preferred name for each of the possible spglib settings """
    output = {}

    with resources.open_text(
            "pyspinw.symmetry.data",
            "preferred_spacegroup_names.txt") as file:
        for line in file:
            parts = line.split(";")
            number = int(parts[0].strip())
            name = parts[1].strip()

            output[number] = name

    return output


formatted_aliases = load_formatted_aliases()
canonical_aliases = load_canonical_aliases()
canonised_to_formatted = load_canonised_to_formatted()
preferred_names = load_preferred_names()