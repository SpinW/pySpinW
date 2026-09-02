import sys
import io
old_stdout = sys.stdout
buffer = io.StringIO()
sys.stdout = buffer
from pyspinw import *
unit_cell = UnitCell(3, 8, 8)
structure = generate_helical_structure(unit_cell, positions=[[0,0,0]], spins=[[0, 1, 0]],
                                   perpendicular=[0,0,1], propagation_vector=[0.5, 0, 0], names=["MCu1"])
exchanges = generate_exchanges(sites=structure,
                               max_distance=3.1,
                               exchange_type=HeisenbergExchange,
                               j=1)
hamiltonian = Hamiltonian(structure, exchanges)
print("#####:0")
hamiltonian.print_summary()
print("######:0")
snapshot(hamiltonian, filename="structure.png")
path = Path([[0,0,0], [1,0,0]])
fig = hamiltonian.spaghetti_plot(path, scale='log', show=False)
fig.savefig("spaghetti_plot.png")
with open("stdout_data.txt", "w") as file:
    file.write(buffer.getvalue())
sys.stdout = old_stdout

