""" Replicates tutorial 20 """
from pyspinw import *

symmetries = '-z, y+3/4, x+3/4; z+3/4, -y, x+3/4; z+3/4, y+3/4, -x; y+3/4, x+3/4, -z; x+3/4, -z, y+3/4; -z, x+3/4, y+3/4'

sg = spacegroup(symmetries) # Gives spacegroup F 4_1/d -3 2/m, second setting, hall=526
print(sg)

a = 10.0307
unit_cell = UnitCell(a, a, a)

yb = LatticeSite(1/2, 1/2, 1/2, 0,0,1/2, name="Yb 3+")
ti = LatticeSite(0,0,0, name="Ti 4+")
o1 = LatticeSite(0.3318, 1/8, 1/8, name="O 2-")
o2 = LatticeSite(3/8, 3/8, 3/8, name="O 2- ")

sites = [yb, ti, o1, o2]

structure = Structure(sites, unit_cell, sg)

# structure.print_summary()
# view(structure)

ybs = structure.sites_by_element("Yb")


# exchanges = generate_exchanges(ybs, unit_cell=unit_cell, bond=1, max_distance=a)
# hamiltonian = Hamiltonian(structure, exchanges)
# hamiltonian.print_summary()
# view(hamiltonian)

structure.print_summary()

base_exchange = HeisenbergExchange(
    structure.site_by_name("Yb 3+"),
    structure.site_by_name("Yb 3+ [1]"),
    j=1,
    name="J")

print("Base")
print(base_exchange)

extra_exchanges = base_exchange.symmetry_fill(structure)

print("By Symmetry")
for exchange in extra_exchanges:
    print(exchange)

symmetry_ham = Hamiltonian(structure, [base_exchange] + extra_exchanges)
view(symmetry_ham)