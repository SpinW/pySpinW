import sys
import io
old_stdout = sys.stdout
old_stderr = sys.stderr
stdout_buffer = io.StringIO()
stderr_buffer = io.StringIO()
sys.stdout = stdout_buffer
sys.stderr = stderr_buffer
from pyspinw import *
unit_cell = UnitCell(1,1,1)
only_site = LatticeSite(1/2, 1/2, 1/2, 0,0,1, name="X")
print("#####:1", file=sys.stderr)
s = Structure([only_site], unit_cell=unit_cell)
print("######:1", file=sys.stderr)
snapshot(s, "ferromagnet_structure.png") #, copies=(5,1,1))
exchanges = [HeisenbergExchange(only_site, only_site, cell_offset=(1,0,0), j=-1)]
exchanges = generate_exchanges(sites=[only_site],
                               unit_cell=unit_cell,
                               max_distance=1.1,
                               exchange_type=HeisenbergExchange,
                               j=-1,
                               direction_filter=filter([1,0,0]))
hamiltonian = Hamiltonian(s, exchanges)
snapshot(s, "ferromagnet_structure.png") #, copies=(5,1,1))
path = Path([[0,0,0], [1,0,0]])
fig = hamiltonian.spaghetti_plot(path, dE=0.4, show=False)
fig.savefig("spaghetti_plot.png")
with open("stdout_data.txt", "w") as file:
    file.write(stdout_buffer.getvalue())
sys.stdout = old_stdout


with open("stderr_data.txt", "w") as file:
    file.write(stderr_buffer.getvalue())
sys.stderr = old_stderr

