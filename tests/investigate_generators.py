from pyspinw.util.group_generators import spglib_generators, parse_one_line_generators
from pyspinw.util.magnetic_symmetry import name_converter
from pyspinw.util.closures import closure

def compare_closure(number):



    generator_string = name_converter.litvin[number].generators

    pyspinw_gen = parse_one_line_generators(generator_string)

    spglib_gen = spglib_generators(number)

    pyspinw_gen.sort()
    spglib_gen.sort()

    remaining_pyspinw = pyspinw_gen.copy()
    pairs = []

    for g1 in spglib_gen:
        for g2 in pyspinw_gen:
            if g1 == g2:
                pairs.append((g1, g2))
                remaining_pyspinw.remove(g2)
                break

        else:
            pairs.append((g1, None))

    for g2 in remaining_pyspinw:
        pairs.append((None, g2))

    # Get string lengths
    left_len = max([len(pair[0].text_form) if pair[0] is not None else 0 for pair in pairs])
    # right_len = max([len(pair[1].text_form) for pair in pairs])

    print(f"{number}:", end="")
    if len(remaining_pyspinw):
        print("Not subset")
    else:
        print("Subset")

    for g1, g2 in pairs:
        s1 = "" if g1 is None else g1.text_form
        s2 = "" if g2 is None else g2.text_form
        s1 += " "*(left_len - len(s1) + 3)

        star = "*" if g1 is not None and g2 is not None else " "
        print(star, s1, s2)




# histogram_comparison()

for i in range(1651):
    compare_closure(i+1)