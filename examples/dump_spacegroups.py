from pyspinw.symmetry.group import database

with open("spacegroup_dump.txt", 'w') as file:
    for i, group in enumerate(database.spacegroups):
        s = f"{i}: {group.number}, {group.choice}, {group.lattice_system}: {group.symbol}\n"
        file.write(s)
