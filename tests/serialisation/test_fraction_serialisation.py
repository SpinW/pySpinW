import pytest
from fractions import Fraction

from pyspinw.serialisation import serialise_fraction, deserialise_fraction, serialise_fraction_or_builtin, \
    deserialise_fraction_or_builtin

fractions = [Fraction(x).limit_denominator() for x in [0.1,0.3,1/3,7/4,18/72,1,7,0]]

@pytest.mark.parametrize("fraction", fractions)
def test_reserialise_fraction(fraction: Fraction):
    serialised = serialise_fraction(fraction)
    deserialised = deserialise_fraction(serialised)

    assert fraction.numerator == deserialised.numerator
    assert fraction.denominator == deserialised.denominator

@pytest.mark.parametrize("fraction", fractions)
def test_reserialise_fraction_or_builtin_fraction(fraction: Fraction):
    serialised = serialise_fraction_or_builtin(fraction)
    deserialised = deserialise_fraction_or_builtin(serialised)

    assert isinstance(deserialised, Fraction)

    assert fraction.numerator == deserialised.numerator
    assert fraction.denominator == deserialised.denominator

@pytest.mark.parametrize("fraction", fractions)
def test_reserialise_fraction_or_builtin_builtin(fraction):
    number = float(fraction)
    serialised = serialise_fraction_or_builtin(number)
    deserialised = deserialise_fraction_or_builtin(serialised)

    assert isinstance(deserialised, float)

    assert number == deserialised

