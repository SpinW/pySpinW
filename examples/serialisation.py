from pyspinw.site import LatticeSite
from pyspinw.coupling import DMCoupling


site1 = LatticeSite(0,0,0, name="A")
site2 = LatticeSite(0.5, 0.5, 0.5, name="B")

coupling = DMCoupling(site1, site2, 1, 1, 1, cell_offset=(0,0,0), name="Example")

json = coupling.serialise()

print(coupling)

print(json)

recovered = coupling.deserialise(json)

print(recovered)
