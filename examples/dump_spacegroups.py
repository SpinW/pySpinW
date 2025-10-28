from pyspinw.symmetry.group import database

with open("spacegroup_dump.txt", 'w') as file:
    for group in database.spacegroups:
        s = f"{group.number}-{group.choice}: {group.symbol}\n"
        file.write(s)