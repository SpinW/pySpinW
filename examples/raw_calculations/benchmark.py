"""Benchmark for each of the raw calculation examples."""

from examples.raw_calculations.ferromagnetic_chain import heisenberg_ferromagnet
from examples.raw_calculations.antiferro_chain import antiferro_chain
from examples.raw_calculations.antiferro_ef import antiferro_ef
from examples.raw_calculations.kagome import kagome_ferromagnet
from examples.raw_calculations.kagome_antiferro import kagome_antiferromagnet
from examples.raw_calculations.kagome_supercell import kagome_supercell
from examples.raw_calculations.utils import InternalClasses, py_classes

from prettytable import PrettyTable
from timeit import timeit

try:
    from pyspinw.rust import spinwave_calculation as rs_spinwave, Coupling as RsCoupling, MagneticField as RsField
    rs_classes = InternalClasses(coupling=RsCoupling, field=RsField)
    RUST_AVAILABLE = True
except ModuleNotFoundError:
    RUST_AVAILABLE = False

from pyspinw.calculations.spinwave import spinwave_calculation as py_spinwave

examples = [
    "heisenberg_ferromagnet",
    "antiferro_chain",
    "antiferro_ef",
    "kagome_ferromagnet",
    "kagome_antiferromagnet",
    "kagome_supercell",
]

rust_option = [False, True] if RUST_AVAILABLE else [False]

for use_rust in rust_option:
    if use_rust:
        spinwave_calculation = "rs_spinwave"
        classes = "rs_classes"
    else:
        spinwave_calculation = "py_spinwave"
        classes = "py_classes"

    table = PrettyTable()
    table.field_names = ["n_q"] + examples

    for n_q in [100, 1000, 5000]:
        times = []
        for example in examples:
            print(f"Benchmarking {example} with n_q={n_q}, rust={use_rust}")
            times.append(
                timeit(
                    f"{spinwave_calculation}(*structure)",
                    setup=f"structure = {example}(n_q={n_q}, classes={classes})",
                    globals=globals(),
                    number=10,
                )
            )

        table.add_row([n_q] + [f"{time:4f}" for time in times])

    print(table)
