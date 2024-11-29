""" Different kinds of anisotropy """

from pyspinw._base import Anisotropy


class GeneralAnisotropy(Anisotropy):
    """ General anisotropy specification """


class DiagonalAnisotropy(Anisotropy):
    """ Anisotropy oriented with axes, but variable amount in x, y and z"""


class XAxisAnisotropy(Anisotropy):
    """ Pure X anisotropy"""


class YAxisAnisotropy(Anisotropy):
    """ Pure Y anisotropy"""


class ZAxisAnisotropy(Anisotropy):
    """ Pure Z anisotropy"""
