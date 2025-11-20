""" This file is responsible for building the table of allowed names

Run this as a script to build the database of name lookups
"""
import os
from collections import defaultdict
from dataclasses import dataclass

from spglib import spglib

from pyspinw.symmetry.settings import Setting, UniqueAxis, RhombohedralOrHexagonal, AxisPermutation


@dataclass
class SettingExpander:
    """ This object is used to convert choices as specified in spglib, where
    defaults are just ommited, to ones where the defaults are explicit"""
    choice_number: bool = False
    permutation: bool = False
    unique_axis: bool = False
    rhombohedral_or_hexagonal: bool = False

    def allow_options_consistent_with(self, setting: Setting):
        """ """
        if setting.choice_number is not None:
            self.choice_number = True

        if setting.permutation is not None:
            self.permutation = True

        if setting.unique_axis is not None:
            self.unique_axis = True

        if setting.rhombohedral_or_hexagonal is not None:
            self.rhombohedral_or_hexagonal = True

    def expand(self, setting: Setting):
        """ Fill in the implicit setting parameters if they are not there """
        return Setting(
            choice_number = 1 if self.choice_number and setting.choice_number is None else setting.choice_number,
            permutation = (AxisPermutation.from_string("abc")
                           if self.permutation and setting.permutation is None
                           else setting.permutation),
            unique_axis = (UniqueAxis.from_string("b")
                           if self.unique_axis and setting.unique_axis is None
                           else setting.unique_axis),
            rhombohedral_or_hexagonal=(RhombohedralOrHexagonal.HEXAGONAL
                                       if self.rhombohedral_or_hexagonal and setting.rhombohedral_or_hexagonal is None
                                       else setting.rhombohedral_or_hexagonal)
        )

def build_expanded_settings() -> dict[int, Setting]:
    """ Create a dictionary that maps to `Setting` objects that have explicit defaults

    If they're part of a group with alternate choices, make sure that there
    is an entry for that kind of choice (e.g. permutation, unique axis) rather than
    None"""
    # Collect together settings for each group
    settings_lists = defaultdict(list)
    for i in range(1, 531):
        group_data = spglib.get_spacegroup_type(i)
        settings_lists[group_data.number].append(i)

    # Set up the objects to control which setting options are applicable to
    # each setting
    expanders: dict[int, SettingExpander] = {} # Map: hall number -> SettingsParameters
    for group_number, settings in settings_lists.items():
        expander = SettingExpander()

        for hall_number in settings:
            setting_data = spglib.get_spacegroup_type(hall_number)
            setting = Setting.from_string(setting_data.choice)
            expander.allow_options_consistent_with(setting)

            expanders[hall_number] = expander

    # Create settings objects for each group
    expanded_settings = {}
    for i in range(1, 531):
        group_data = spglib.get_spacegroup_type(i)
        expanded_settings[i] = expanders[i].expand(Setting.from_string(group_data.choice))

    return expanded_settings

def build_spacegroup_lookup(target_directory: str):
    """
    Create lookup tables for getting spacegroups from strings

    See ADR 010 for requirements.
    """

    default_names = {}
    default_setting_index = {}
    full_names = {}
    bracket_names = {}

    expanded_settings = build_expanded_settings()

    for i in range(1, 531):
        group_data = spglib.get_spacegroup_type(i)

        setting = Setting.from_string(group_data.choice)
        expanded_setting = expanded_settings[i]

        # First group by number should be the default
        if group_data.number not in default_setting_index:
            default_names[group_data.number] = group_data.international_short
            default_setting_index[group_data.number] = i

        has_non_default = False
        # Add full name, for the first entry only
        if group_data.international_full not in full_names:
            full_names[group_data.international_full] = i
            has_non_default = True

        # Groups 67 and 68 need some special treatment because of permuations
        if group_data.number == 67 or group_data.number == 68:
            if setting.permutation is not None:
                bracket_names[group_data.international_full + f" ({setting.permutation.to_string()})"] = i
                bracket_names[group_data.international_full + f" [{setting.permutation.to_string()}]"] = i
                has_non_default = True

        # Add choice number if choice number is the only variation
        if (expanded_setting.choice_number is not None and
            expanded_setting.unique_axis is None and
            expanded_setting.permutation is None and
            expanded_setting.rhombohedral_or_hexagonal is None):

            full_names[group_data.international_full + " : "+str(setting.choice_number)] = i
            has_non_default = True

        # Add names with R or H suffix
        if (expanded_setting.choice_number is None and
            expanded_setting.unique_axis is None and
            expanded_setting.permutation is None and
            expanded_setting.rhombohedral_or_hexagonal is not None):

            r_or_h = " " + expanded_setting.rhombohedral_or_hexagonal.value

            full_names[group_data.international_full + r_or_h ] = i
            full_names[group_data.international_short + r_or_h] = i

            has_non_default = True

        # Add short name with choice in brackets
        if group_data.choice != "":
            setting_string = expanded_setting.to_string()
            prefix = default_names[group_data.number]
            bracket_names[prefix + f" ({setting_string})"] = i
            bracket_names[prefix + f" [{setting_string}]"] = i
            has_non_default = True



    if not has_non_default:
        raise Exception("Setting", i, "has no non-default entry")
    #
    # print(len(default_setting_index), "defaults")

    # Build list of names for each group, we want it so
    # that the first entry is the one we want to use to identify the group
    full_lookup = defaultdict(list)
    for group_number in default_setting_index:
        full_lookup[default_setting_index[group_number]].append(default_names[group_number])

    for name, number in full_names.items():
        full_lookup[number].append(name)

    for name, number in bracket_names.items():
        full_lookup[number].append(name)

    for i in range(1, 531):
        if i not in full_lookup:
            print(f"No entry for setting {i}")

    preferred_names = defaultdict(str)
    for number, name_list in full_lookup.items():
        if not name_list:
            raise Exception("This should never be reached, because there should not be any "
                            "entry in the dict if the list is empty")
        preferred_names[number] = name_list[0]

    # Check that the preferred names are unique

    if len(set([key for key in preferred_names.keys()])) != 530:
        raise Exception("Preferred name collision (two or more groups/settings have the"
                        "same name)")

    # Write out the preferred names
    with open(os.path.join(target_directory, "preferred_spacegroup_names.txt"), 'w') as file:
        for i in range(1, 531):
            file.write(f"{i-1}: ")
            file.write(preferred_names[i])
            file.write("\n")


    with open(os.path.join(target_directory, "all_spacegroup_names.txt"), 'w') as file:
        for i in range(1, 531):
            group_data = spglib.get_spacegroup_type(i)
            file.write(f"{i}; {group_data.number}; ")
            file.write("; ".join(full_lookup[i]))
            file.write("\n")


if __name__ == "__main__":

    build_spacegroup_lookup(target_directory="data")
