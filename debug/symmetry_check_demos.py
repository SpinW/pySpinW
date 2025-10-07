from pyspinw.oldstyle.supercell import rotation_supercell
from pyspinw.site import LatticeSite
from pyspinw.symmetry.group import magnetic_group_symbol_lookup
from pyspinw.symmetry.supercell import TrivialSupercell, TransformationSupercell, CommensuratePropagationVector, \
    RotationTransform
from pyspinw.symmetry.symmetry_checking import check_supercell_moment_consistency

def display_warnings(warnings):
    for offset in warnings:
        print(f"Cell {offset}: {len(warnings[offset])} warnings")

    for offset in warnings:
        for warning in warnings[offset]:
            print(warning)


def check_moment_consistency_p1_example():
    sites = [LatticeSite(0.5,0.5,0.5, 0, 0, 1)]
    symmetry = magnetic_group_symbol_lookup["P1.1"]
    supercell = TransformationSupercell([
        (CommensuratePropagationVector(1/3, 0, 0), RotationTransform([0, 1, 0])),
        (CommensuratePropagationVector(0, 1/3, 0), RotationTransform([0, 1, 0])),
        (CommensuratePropagationVector(0, 0, 1/3), RotationTransform([0, 1, 0]))])

    warnings = check_supercell_moment_consistency(supercell, symmetry, sites)

    display_warnings(warnings)
    print("There should not have been any warnings listed")


def check_moment_consistency_forbidden_example(group="P4_2/mnm.1'_I[I4/mmm]"):
    sites = [LatticeSite(0.5,0.5,0.5, 0, 0, 1),
             LatticeSite(0.1,0.1,0.1,0,0,1)]
    symmetry = magnetic_group_symbol_lookup[group]
    supercell = TransformationSupercell([
        (CommensuratePropagationVector(1/3, 0, 0), RotationTransform([0, 1, 0])),
        (CommensuratePropagationVector(0, 1/3, 0), RotationTransform([1, 0, 0])),
        (CommensuratePropagationVector(0, 0, 1/3), RotationTransform([0, 1, 0]))])

    warnings = check_supercell_moment_consistency(supercell, symmetry, sites)

    display_warnings(warnings)


def find_supercell_only_inconsistent_groups():
    """ List all the groups that are inconsistent with the supercell, but not the base unit cell,
    for one particular system"""

    sites = [LatticeSite(0,0,0,1,0,0)]
    supercell = TransformationSupercell([
        (CommensuratePropagationVector(1/3, 0, 0), RotationTransform([0, 1, 0])),
        (CommensuratePropagationVector(0, 1/3, 0), RotationTransform([1, 0, 0])),
        (CommensuratePropagationVector(0, 0, 1/3), RotationTransform([0, 1, 0]))])

    for group_name in magnetic_group_symbol_lookup:
        # print("Checking",group_name,"...")

        group = magnetic_group_symbol_lookup[group_name]
        warnings = check_supercell_moment_consistency(supercell, group, sites)

        n = 0
        n_0 = 0
        for offset in warnings:
            if offset == (0,0,0):
                n_0 = len(warnings[offset])
            else:
                n = max([len(warnings[offset]), n])

        if n_0 == 0 and n != 0:
            print(group_name)
            # display_warnings(warnings)


def check_moment_consistency_forbidden_example_2(group="P4_2/mnm.1'_I[I4/mmm]"):
    sites = [LatticeSite(0.5,0.5,0.5, 0, 0, 1)]
    symmetry = magnetic_group_symbol_lookup[group]
    supercell = TransformationSupercell([
        (CommensuratePropagationVector(1/3, 0, 0), RotationTransform([0, 1, 0])),
        (CommensuratePropagationVector(0, 1/3, 0), RotationTransform([1, 0, 0])),
        (CommensuratePropagationVector(0, 0, 1/3), RotationTransform([0, 1, 0]))])

    warnings = check_supercell_moment_consistency(supercell, symmetry, sites)


    for offset in warnings:
        for warning in warnings[offset]:
            print(f"Cell {offset}: ", end='')
            print(warning)

if __name__ == "__main__":

    check_moment_consistency_p1_example()

    # check_moment_consistency_forbidden_example()

    # find_supercell_only_inconsistent_groups()
