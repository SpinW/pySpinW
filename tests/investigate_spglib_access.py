import pytest
import spglib
import numpy as np

from pyspinw.util.magnetic_symmetry import name_converter
from pyspinw.util.group_generators import spglib_generators

def give_all_symmetries(number):
    print("BNS:", name_converter.litvin[number].bns_symbol)
    print("UNI:", name_converter.litvin[number].uni_symbol)

    gens = spglib_generators(number)

    for gen in gens:
        print(gen.text_form[1:-1])

give_all_symmetries(65)