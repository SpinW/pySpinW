"""Jupyter Notebook visualization for PySpinW using ipywidgets and ipywebgl."""

from typing import Tuple

_REQUIRED_PACKAGES = {
    "ipywidgets": "ipywidgets",
    "ipywebgl": "ipywebgl",
    "ipyevents": "ipyevents",
}


def check_jupyter_dependencies() -> None:
    """Verify that all optional dependencies required for Jupyter visualization are installed.

    Raises
    ------
    ImportError
        If one or more of ``ipywidgets``, ``ipywebgl``, or ``ipyevents`` is missing.
    """
    missing = []
    for mod_name, pkg_name in _REQUIRED_PACKAGES.items():
        try:
            __import__(mod_name)
        except ImportError:
            missing.append(pkg_name)

    if missing:
        missing_str = ", ".join(missing)
        raise ImportError(
            f"Missing optional dependencies for Jupyter visualization: {missing_str}.\n"
            f"Please install them using: pip install 'pyspinw[jupyter]' or pip install {' '.join(missing)}"
        )


from pyspinw.gui.jupyter.interaction import CanvasInteractionHandler
from pyspinw.gui.jupyter.crystalview_jupyter import CrystalViewerWidget


__all__ = ["check_jupyter_dependencies", "CanvasInteractionHandler", "CrystalViewerWidget"]

