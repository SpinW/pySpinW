import sys
import io
old_stdout = sys.stdout
old_stderr = sys.stderr
stdout_buffer = io.StringIO()
stderr_buffer = io.StringIO()
sys.stdout = stdout_buffer
sys.stderr = stderr_buffer


from pyspinw import *



two_site_unit_cell = UnitCell(2, 1, 1)
two_site_site_1 = LatticeSite(0.25, 0.5, 0.5, sz=1, name="a")
two_site_site_2 = LatticeSite(0.75, 0.5, 0.5, sz=-1, name="b")
two_site_structure = Structure([two_site_site_1, two_site_site_2], two_site_unit_cell)

snapshot(two_site_structure, "two_site_structure.png", view_point=(0,-5,-0.5))


two_site_exchanges = [
    HeisenbergExchange(two_site_site_1, two_site_site_2, j=1),
    HeisenbergExchange(two_site_site_2, two_site_site_1, cell_offset=(1,0,0), j=1)
]

two_site_exchanges = generate_exchanges([two_site_site_1, two_site_site_2],
                               unit_cell=two_site_unit_cell,
                               max_distance=1.1,
                               direction_filter=filter([1,0,0], symmetric=True))

print("#####:1")
for exchange in two_site_exchanges:
    print(exchange)
print("######:1")

two_site_hamiltonian = Hamiltonian(two_site_structure, two_site_exchanges)
path = Path([(0, 0, 0), (1, 0, 0)], convert_to_lattice_units_with=two_site_unit_cell)
fig = two_site_hamiltonian.spaghetti_plot(path, show=False)
fig.savefig("two_site_dispersion.png")

unit_cell = UnitCell(1,1,1)

structure = generate_helical_structure(unit_cell, positions=[[0,0,0]], spins=[[0, 1, 0]],
                                   perpendicular=[0,0,1], propagation_vector=[0.5, 0, 0], names=["MCu1"])


exchanges = generate_exchanges(sites=structure,
                               max_distance=3.1,
                               exchange_type=HeisenbergExchange,
                               j=1)

hamiltonian = Hamiltonian(structure, exchanges)

path = Path([(0,0,0), (1,0,0)])

print("#####:2")
hamiltonian.print_summary()
print("######:2")


snapshot(hamiltonian, filename="structure.png", view_point=(-5,-5,-10),
         display_options=DisplayOptions(perspective=False, atom_spin_scaling=0.5))

path = Path([[0,0,0], [1,0,0]])


fig = hamiltonian.spaghetti_plot(path, scale='log', show=False)
fig.savefig("spaghetti_plot.png")

with open("stdout_data.txt", "w") as file:
    file.write(stdout_buffer.getvalue())
sys.stdout = old_stdout


with open("stderr_data.txt", "w") as file:
    file.write(stderr_buffer.getvalue())
sys.stderr = old_stderr

