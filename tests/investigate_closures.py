from pyspinw.util.group_generators import spglib_generators, parse_one_line_generators
from pyspinw.util.magnetic_symmetry import name_converter
from pyspinw.util.closures import closure

def compare_closure(number):



    generator_string = name_converter.litvin[number].generators

    pyspinw_gen = closure(parse_one_line_generators(generator_string))

    spglib_gen = spglib_generators(number)

    pyspinw_gen.sort()
    spglib_gen.sort()

    print(f"{number}: ", end="")

    if len(pyspinw_gen) != len(spglib_gen):
        print("Different group size")
        return

    mismatches = []
    for g1, g2 in zip(pyspinw_gen, spglib_gen):
        if g1 != g2:
            mismatches.append((g1, g2))

    if mismatches:
        print("Mismatch")
        for g1, g2 in mismatches:
            print(g1, "vs", g2)

        return

    else:
        print("Same")

for i in range(1651):
    compare_closure(i+1)