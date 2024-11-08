""" Different Kinds of samples"""

from pyspinw._base import Sample


class SingleCrystal(Sample):
    """ Specifies a single crystal sample """


class Multidomain(Sample):
    """ Sample consisting of multiple domains"""


class Twin(Multidomain):
    """ Specify a twinned crystal.

    Special case of multidomain """


class Powder(Sample):
    """ Sample is a powder """


