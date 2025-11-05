"""Utility functions for the raw calculation examples."""

from dataclasses import dataclass
from typing import Callable

from pyspinw.calculations.spinwave import (
    spinwave_calculation as py_spinwave,
    Coupling as PyCoupling,
    MagneticField as PyField,
)

@dataclass
class InternalClasses:
    """A system of internal classes."""

    coupling: object
    field: object

py_classes = InternalClasses(coupling=PyCoupling, field=PyField)


def run_example(
    example: Callable[[int, object], tuple[list, "np.array", "np.array", list]],
    use_rust: bool = False,
) -> (tuple, "np.array"):
    """Run an example.

    Parameters
    ----------
    example: Callable
        A function which takes two arguments:

        - n_q: the number of q-values to calculate for
        - coupling_class: the class to use for coupling objects

        and returns four objects:

        - a list of rotations;
        - an array of magnitudes;
        - an array of q-vectors;
        - a list of couplings.

    use_rust: bool
        Whether to use the Rust calculation routine.

    Returns
    -------
    structure: tuple
        The structure from the example callable.
    energies: np.array
        The result from the spinwave calculation.

    """
    try:
        from pyspinw.rust import spinwave_calculation as rs_spinwave, Coupling as RsCoupling, MagneticField as RsField

        RUST_AVAILABLE = True
    except ModuleNotFoundError:
        RUST_AVAILABLE = False

    # default to Python unless Rust is requested (which it is by default) and available
    classes = py_classes
    spinwave_calculation = py_spinwave
    if use_rust:
        if RUST_AVAILABLE:
            classes = InternalClasses(coupling=RsCoupling, field=RsField)
            spinwave_calculation = rs_spinwave
        else:
            print(
                "Attempted to run example in Rust, but Rust is not available. "
                "To run in Python and hide this warning, "
                "pass 'python' as an argument to the Python script in the terminal."
            )

    structure = example(classes=classes)
    energies, sab = spinwave_calculation(*structure)

    return structure, energies, sab
