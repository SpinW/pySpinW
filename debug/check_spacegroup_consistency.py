from pyspinw.symmetry.group import database


def check_consistency(spacegroup):

    for op1 in spacegroup.operations:
        for op2 in spacegroup.operations:
            found_any = False
            for op3 in spacegroup.operations:
                found_any = found_any or op1.and_then(op2) == op3

            assert found_any, f"{spacegroup}: Composition of {op1} and {op2} not in group {op1.and_then(op2)}"

if __name__ == "__main__":



    with open("spacegroup_checks.txt", 'w') as file:

        def printy(s):
            print(s)
            file.write(f"{s}\n")

        for spacegroup in database.spacegroups:
            try:
                check_consistency(spacegroup)
                printy(f"{spacegroup}: OK")
            except AssertionError as ae:
                printy(ae)
