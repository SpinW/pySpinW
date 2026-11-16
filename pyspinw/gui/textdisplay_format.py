"""Text format utilities for the TextDisplay class without Qt/ipy imports"""

import numpy as np
from pyspinw.anisotropy import Anisotropy

def format_triple(triple):
    """Rendering for numeric triples.

    Parameters
    ----------
    triple : tuple or array_like of float
        Sequence of three numeric values.

    Returns
    -------
    str
        Formatted string representation ``"(a, b, c)"``.
    """
    a,b,c = triple
    return f"({a:.3g}, {b:.3g}, {c:.3g})"

def format_spin_data(spin_data):
    """Rendering for spin data.

    Parameters
    ----------
    spin_data : numpy.ndarray
        2D array of shape ``(n, 3)`` containing spin vectors.

    Returns
    -------
    str
        Comma-separated string of formatted spin triples.
    """
    return ", ".join([format_triple(spin_data[i, :]) for i in range(spin_data.shape[0])])

def format_anisotropies(anisotropies: list[Anisotropy]):
    """Formatting for the anisotropies column.

    Parameters
    ----------
    anisotropies : list of Anisotropy
        List of anisotropy terms associated with a site.

    Returns
    -------
    str
        Semicolon-separated parameter strings, or ``"None"`` if empty.
    """
    if anisotropies:
        return "; ".join([a.parameter_string for a in anisotropies])
    else:
        return "None"

def format_g(g_matrix: np.ndarray):
    """Formatting for g-factors.

    Parameters
    ----------
    g_matrix : (3, 3) numpy.ndarray
        3x3 g-tensor matrix.

    Returns
    -------
    str
        Compact string representation of the g-tensor.
    """
    diag = np.diagonal(g_matrix)
    if np.all(g_matrix == np.diag(diag)):
        a,b,c = diag

        if a==b and b==c:
            return f"{a:.3g}"
        else:
            return f"{a:.3g}, {b:.3g}, {c:.3g}"
    else:
        m = g_matrix.reshape(-1)
        return (f"[[{m[0]:.3g}, {m[1]:.3g}, {m[2]:.3g}],"
                f" [{m[3]:.3g}, {m[4]:.3g}, {m[5]:.3g}],"
                f" [{m[6]:.3g}, {m[7]:.3g}, {m[8]:.3g}]]")
