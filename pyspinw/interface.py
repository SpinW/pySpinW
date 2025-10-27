""" Helper functions for python interface """

from pyspinw.symmetry.group import database

def spacegroup(name: str):
    return database.spacegroup_by_name(name)