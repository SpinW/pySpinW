from pyspinw import *

p4mmm = spacegroup("p4/mmm")

centre = LatticeSite(0.5, 0.5, 0.5, name="centre")
offset = LatticeSite(0.5, 0.5, 0.1, name="offset")

cell = UnitCell(1,1,1)

struct = Structure(sites=[centre, offset], unit_cell=cell, spacegroup=p4mmm)

struct.print_summary()

other_offset = struct.site_by_name("offset [1]")
exchange_1 = DMExchange(centre, offset, d_x=1, d_y=1, d_z=1)

exchange_2 = exchange_1.symmetry_copy(centre, other_offset, p4mmm)

print(exchange_2)