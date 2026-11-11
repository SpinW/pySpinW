"""Visualisation and crystal structure editing GUI utilities."""

from pyspinw.hamiltonian import Hamiltonian
from pyspinw.structure import Structure

try:
    from pyspinw.gui.viewer import show_object, snapshot
    HASQT = True
except ModuleNotFoundError:
    HASQT = False
    def snapshot(*args, **kwargs):
        raise ModuleNotFoundError('Could not import PySide6/OpenGL for snapshot')

try:
    from pyspinw.gui.jupyter.viewer import JupyterViewer
    HASJUPYTER = True
except ModuleNotFoundError:
    HASJUPYTER = False


def view(object: Hamiltonian | Structure,
         block=True,
         copies: tuple[int,int,int]=(1,1,1),
         rotation_supercell_expansion_max: tuple[int,int,int]=(10,10,10),
         propagation_vector_commensurability_tolerance=1e-4,
         jupyter: bool=False):
    """Show a Hamiltonian or structure in a viewer within the Jupyter notebook

    Parameters
    ----------
    object : Hamiltonian or Structure
        Hamiltonian or Structure object to display.
    block : bool, optional
        Whether to block execution by running the Qt event loop. Default is
        ``True``.
    copies : tuple of (int, int, int), optional
        Supercell scaling along each lattice vector. Default is ``(1, 1, 1)``.
    rotation_supercell_expansion_max : tuple of (int, int, int), optional
        Maximum denominator along each axis when approximating a
        ``RotationSupercell`` propagation vector as a commensurate cell.
        Default is ``(10, 10, 10)``.
    propagation_vector_commensurability_tolerance : float, optional
        Tolerance for matching a propagation vector component to a rational
        fraction. Default is ``1e-4``.
    jupyter : bool, optional
        Whether to display Hamiltonian or Structure in a Jupyter widget
        or use the default external window. Default is ``False``.

    Raises
    ------
    TypeError
        If ``object`` is not a ``Hamiltonian`` or ``Structure``.
    """
    if jupyter:
        if HASJUPYTER:
            display(JupyterViewer(object, copies, rotation_supercell_expansion_max,
                propagation_vector_commensurability_tolerance))
            return
        else:
            raise ModuleNotFoundError('Could not import ipywebgl or ipywidgets libraries')
    if not HASQT:
        raise ModuleNotFoundError('Could not import PySide or OpenGL libraries')

    show_object(object, block, copies, rotation_supercell_expansion_max,
                propagation_vector_commensurability_tolerance)
