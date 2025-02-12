from Tools.scripts.stable_abi import generators

from pyspinw.util.group_generators import spglib_generators, parse_one_line_generators
from pyspinw.util.magnetic_symmetry import name_converter


def check_subset(number):
    print(f"{number}:", end="")
    group = spglib_generators(number)
    gens = parse_one_line_generators(name_converter.litvin[number].generators)

    all_found = True
    misses = []
    for generator in gens:
        found = False
        for symmetry in group:
            if symmetry == generator:
                found = True
                break

        if not found:
            misses.append(generator)

        all_found = all_found and found

    print("Subset" if all_found else "")

    for miss in misses:
        print("  ", miss.text_form)


def print_generators(number):
    print(name_converter.litvin[number].bns_symbol)

    gens = parse_one_line_generators(name_converter.litvin[number].generators)

    for gen in gens:
        print(gen.text_form)


for i in range(1, 1652):
    check_subset(i)
#
# print_generators(65)