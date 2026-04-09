""" Spin wave spectrum of Ba3NbFe3Si2O14 """

from pyspinw import *

unit_cell = UnitCell(8.539,8.539, 5.2414, 90, 90, 120)
sites = [LatticeSite(0.24964,0, 1/2, 0, 0, 5/2)]

structure = Structure(sites, unit_cell, spacegroup("P321"))



view(structure)

#exchanges = generate_exchanges(sites)