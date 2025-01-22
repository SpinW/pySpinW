import pytest
import spglib

from pyspinw.util.magnetic_symmetry import name_converter
from pyspinw.util.group_generators import parse_one_line_generators



@pytest.mark.parametrize("number", range(1, 1652))
def test_generator_consistency(number: int):
    generator_string = name_converter.litvin[number].generators

    generators_from_pyspinw_database = parse_one_line_generators(generator_string)


    generators_from_spglib = spglib.get_symmetry_from_database(number)

    print(generators_from_spglib)
    print(generators_from_pyspinw_database)