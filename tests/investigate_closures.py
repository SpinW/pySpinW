from pyspinw.util.group_generators import spglib_generators, parse_one_line_generators
from pyspinw.util.magnetic_symmetry import name_converter
from pyspinw.util.closures import closure

def compare_closure(number):



    generator_string = name_converter.litvin[number].generators

    initial_generators = parse_one_line_generators(generator_string)
    pyspinw_gen = closure(initial_generators)

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

    print(f"{number} ({name_converter.litvin[number].bns_number}):", end="")
    if any([pair[0] is None or pair[1] is None for pair in pairs]):
        print("Not matching")

        name1 = "spglib"
        name2 = "spinw"

        print(" ", (" "*(left_len - len(name1))) + name1, " ", name2)

        for g1, g2 in pairs:
            s1 = "" if g1 is None else g1.text_form
            s2 = "" if g2 is None else g2.text_form
            s1 = " "*(left_len - len(s1)) + s1

            equal = "=" if g1 is not None and g2 is not None else " "

            in_generators = False
            if g2 is not None:
                for initial in initial_generators:
                    in_generators = in_generators or initial == g2

            star = "*" if in_generators else ""
            print(" ", s1, equal, s2, star)


    else:
        print("Same")


def histogram_comparison():

    spw = []
    spg = []

    for number in range(1, 1652):
        print(number)

        generator_string = name_converter.litvin[number].generators

        pyspinw_gen = closure(parse_one_line_generators(generator_string))
        spw.append(len(pyspinw_gen))

        spglib_gen = spglib_generators(number)
        spg.append(len(spglib_gen))

    import matplotlib.pyplot as plt

    plt.scatter(spw, spg)
    plt.show()



# histogram_comparison()

for i in range(1651):
    compare_closure(i+1)