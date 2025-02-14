from pyspinw.util.closures import closure
from pyspinw.util.group_generators import spglib_generators, parse_one_line_generators
from pyspinw.util.magnetic_symmetry import name_converter

print("Loading data")
spglib_operators = [spglib_generators(i) for i in range(1, 1652)]
spinw_generators = [parse_one_line_generators(name_converter.litvin[i].generators) for i in range(1, 1652)]

def match_spinw(number, spinw_gens):
    matching = []
    spinw_ops = closure(spinw_gens)
    for spglib_number, spglib_ops in enumerate(spglib_operators):
        spglib_ops.sort()
        spinw_ops.sort()

        if len(spinw_ops) != len(spglib_ops):
            continue

        match = True
        for a, b in zip(spinw_ops, spglib_ops):
            if a != b:
                match = False
                break

        if match:
            matching.append(spglib_number+1)


    if len(matching) == 0:
        print(number+1, "no match,", name_converter.litvin[number].bns_symbol, name_converter.litvin[number].uni_symbol)
        print("  generators:")
        for spinw_gen in spinw_gens:
            print("    ", spinw_gen)
    else:
        print(number+1, "matches", matching)

for i in range(1651):
    match_spinw(i, spinw_generators[i])

