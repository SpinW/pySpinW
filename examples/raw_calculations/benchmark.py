"""Benchmark for each of the raw calculation examples."""
from examples.raw_calculations.ferromagnetic_chain import heisenberg_ferromagnet
from examples.raw_calculations.antiferro_chain import antiferro_chain
from examples.raw_calculations.kagome import kagome_ferromagnet
from examples.raw_calculations.kagome_supercell import kagome_supercell

from prettytable import PrettyTable
from timeit import timeit

from pyspinw.calculations.spinwave import spinwave_calculation

table = PrettyTable()
table.field_names = ["n_q", "ferromagnetic_chain", "antiferro_chain", "kagome_ferromagnet", "kagome_supercell"]

for n_q in [100, 1000, 5000]:
    times = []
    for example in ["heisenberg_ferromagnet", "antiferro_chain", "kagome_ferromagnet", "kagome_supercell"]:
        times.append(timeit('spinwave_calculation(*structure)', setup=f'structure = {example}(n_q={n_q})', globals=globals(), number=10))

    table.add_row([n_q] + [f"{time:4f}" for time in times])

print(table)

