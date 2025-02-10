import pytest
import spglib
import numpy as np

from pyspinw.util.magnetic_symmetry import name_converter
from pyspinw.util.group_generators import parse_one_line_generators, _spglib_generators_to_objects
from pyspinw.util.apply_generators import apply_generators_with_moments, apply_generators_until_stable

# np.random.seed(1234)
test_points = np.random.random((4, 6))

def generator_consistency(number: int):
    generator_string = name_converter.litvin[number].generators

    generators_from_pyspinw_database = parse_one_line_generators(generator_string)

    generators_from_spglib = _spglib_generators_to_objects(spglib.get_magnetic_symmetry_from_database(number))

    #
    # print(generators_from_spglib)
    # print(generators_from_pyspinw_database)

    n_spglib = len(generators_from_spglib)
    # n_other = len(generators_from_spglib["rotations"]), len(generators_from_spglib["translations"])


    print(f"{number}: {len(generators_from_pyspinw_database)}, {n_spglib}")

    return generators_from_spglib, generators_from_pyspinw_database

def numerical_inspect(number: int):
    generator_string = name_converter.litvin[number].generators

    generators_from_pyspinw_database = parse_one_line_generators(generator_string)


    generators_from_spglib = _spglib_generators_to_objects(spglib.get_magnetic_symmetry_from_database(number))

    spinw_points = apply_generators_until_stable(test_points, generators_from_pyspinw_database)
    spglib_points = apply_generators_with_moments(test_points, generators_from_spglib)

    print(spinw_points.shape, len(generators_from_pyspinw_database), spglib_points.shape, len(generators_from_spglib))

def full_sort(points):
    for i in range(6):
        points = points[np.argsort(points[:, i], kind="mergesort"), :]
    return points


def numerical_check(number: int, tol=1e-8):
    generator_string = name_converter.litvin[number].generators

    generators_from_pyspinw_database = parse_one_line_generators(generator_string)


    generators_from_spglib = _spglib_generators_to_objects(spglib.get_magnetic_symmetry_from_database(number))

    # print("Applying generators from SpinW database")
    spinw_points = apply_generators_until_stable(test_points, generators_from_pyspinw_database)

    # print("Applying generators from spglib")
    spglib_points = apply_generators_with_moments(test_points, generators_from_spglib)

    spinw_points = full_sort(spinw_points)
    spglib_points = full_sort(spglib_points)


    if spinw_points.shape != spglib_points.shape:
        print(f"{number}: Different ({spinw_points.shape[0]} vs {spglib_points.shape[0]})")
        return "count error"
        # print(spinw_points)
        # print(spglib_points)
    else:
        n = spinw_points.shape[0] * spinw_points.shape[1]
        delta_sq = np.sum((spinw_points - spglib_points) ** 2) / n
        if  delta_sq > tol * tol:
            # import matplotlib.pyplot as plt
            # plt.subplot(1,2,1)
            # plt.scatter(spinw_points[:, 0], spinw_points[:, 1])
            # plt.scatter(spglib_points[:, 0], spglib_points[:, 1])
            #
            # plt.subplot(1,2,2)
            # plt.scatter(spinw_points[:, 0], spinw_points[:, 2])
            # plt.scatter(spglib_points[:, 0], spglib_points[:, 2])
            #
            # plt.show()
            #
            # print("\n\n\nSpinW")
            # print(spinw_points)
            # print("\n\n\nspglib")
            # print(spglib_points)

            print(f"{number}: Different (delta_sq={delta_sq})")
            return "numerical error"
        else:
            print(f"{number}: Same")
            return "same"




# for i in range(1, 1652):
#     generator_consistency(i)
#
#
#
# spglib_generators, pyspinw_generators = generator_consistency(3)
#
# print(pyspinw_generators)
# print(spglib_generators)
#
# n = 130
# print(name_converter.litvin[n].generators)
# print(name_converter.litvin[n].bns_number)
# numerical_check(n)

with open("matching_groups.txt", 'w') as match:
    with open("non_matching_groups_numerical.txt", 'w') as nonmatch_numerical:
        with open("non_matching_groups_count.txt", 'w') as nonmatch_count:

            for i in range(1, 1652):

                match numerical_check(i):
                    case "numerical error":
                        writer = nonmatch_numerical
                    case "count error":
                        writer = nonmatch_count
                    case "same":
                        writer = match
                    case _:
                        raise Exception("You done gone wrong")

                writer.write(f"{name_converter.litvin[i].bns_symbol}\n")

# numerical_check(5)
# numerical_check(1641)