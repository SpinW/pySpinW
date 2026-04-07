from pyspinw.site import LatticeSite
from pyspinw.exchange import DMExchange


site1 = LatticeSite(0,0,0, name="A")
site2 = LatticeSite(0.5, 0.5, 0.5, name="B")

exchange = DMExchange(site1, site2, 1, 1, 1, cell_offset=(0, 0, 0), name="Example")

json = exchange.serialise()

print(exchange)

print(json)

recovered = exchange.deserialise(json)

print(recovered)
