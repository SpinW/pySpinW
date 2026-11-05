"""Composite Jupyter Notebook viewer combining WebGL canvas, toolbar controls, and model metadata."""

from __future__ import annotations

import ipywidgets as widgets
import numpy as np

from pyspinw.gui.displayoptions import DisplayOptions
from pyspinw.gui.jupyter import check_jupyter_dependencies
from pyspinw.gui.jupyter.controls import DisplayControlsWidget
from pyspinw.gui.jupyter.crystalview_jupyter import CrystalViewerWidget
from pyspinw.gui.rendermodel import RenderModel
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.structure import Structure


class ViewerJupyter(widgets.VBox):
    """Composite Jupyter Notebook widget combining WebGL 3D viewer and controls.

    Parameters
    ----------
    hamiltonian : Hamiltonian
        The magnetic Hamiltonian or structure to display.
    initial_rotation : (3, 3) numpy.ndarray or None, optional
        Initial 3x3 view rotation matrix. Default is None.
    initial_distance : float or None, optional
        Initial camera distance from the view origin. Default is None.
    display_options : DisplayOptions or None, optional
        Display options. Default is None.
    copies : tuple of (int, int, int), optional
        Supercell scaling along each lattice vector. Default is (1, 1, 1).
    rotation_supercell_expansion_max : tuple of (int, int, int), optional
        Maximum denominator along each axis when approximating a
        RotationSupercell propagation vector as a commensurate cell. Default is (10, 10, 10).
    propagation_vector_commensurability_tolerance : float, optional
        Tolerance for matching a propagation vector component to a rational fraction. Default is 1e-4.
    width : int, optional
        Canvas width in pixels. Default is 700.
    height : int, optional
        Canvas height in pixels. Default is 500.
    """

    def __init__(
        self,
        hamiltonian: Hamiltonian,
        initial_rotation: np.ndarray | None = None,
        initial_distance: float | None = None,
        display_options: DisplayOptions | None = None,
        copies: tuple[int, int, int] = (1, 1, 1),
        rotation_supercell_expansion_max: tuple[int, int, int] = (10, 10, 10),
        propagation_vector_commensurability_tolerance: float = 1e-4,
        width: int = 700,
        height: int = 500,
    ):
        check_jupyter_dependencies()
        super().__init__()

        self.hamiltonian = hamiltonian
        self.render_model = RenderModel(
            hamiltonian=hamiltonian,
            copies=copies,
            rotation_supercell_expansion_max=rotation_supercell_expansion_max,
            propagation_vector_tolerance=propagation_vector_commensurability_tolerance,
        )

        # 1. 3D WebGL Canvas
        self.viewer = CrystalViewerWidget(
            render_model=self.render_model,
            initial_rotation=initial_rotation,
            initial_distance=initial_distance,
            display_options=display_options,
            width=width,
            height=height,
        )

        # 2. Controls Toolbar
        self.controls = DisplayControlsWidget(viewer=self.viewer)

        # 3. Info / Status banner
        uc = self.hamiltonian.structure.unit_cell
        n_sites = len(self.hamiltonian.structure.sites)
        n_mag = sum(1 for s in self.render_model.sites if s.is_magnetic)
        n_bonds = len(self.hamiltonian.exchanges)

        info_text = (
            f"<b>Lattice:</b> a={uc.a:.3f}Å, b={uc.b:.3f}Å, c={uc.c:.3f}Å, "
            f"α={uc.alpha:.1f}°, β={uc.beta:.1f}°, γ={uc.gamma:.1f}° | "
            f"<b>Sites:</b> {n_sites} ({n_mag} magnetic) | "
            f"<b>Bonds:</b> {n_bonds}"
        )
        self.info_label = widgets.HTML(
            value=f"<div style='font-size: 12px; color: #444; padding: 4px 2px; font-family: sans-serif;'>{info_text}</div>"
        )

        # Assemble layout
        self.children = [self.info_label, self.controls, self.viewer]
        self.layout.border = "1px solid #ddd"
        self.layout.padding = "8px"
        self.layout.border_radius = "4px"


def view_jupyter(
    target: Hamiltonian | Structure,
    initial_rotation: np.ndarray | None = None,
    initial_distance: float | None = None,
    display_options: DisplayOptions | None = None,
    copies: tuple[int, int, int] = (1, 1, 1),
    rotation_supercell_expansion_max: tuple[int, int, int] = (10, 10, 10),
    propagation_vector_commensurability_tolerance: float = 1e-4,
    width: int = 700,
    height: int = 500,
) -> ViewerJupyter:
    """Display a crystal structure or Hamiltonian in an interactive 3D Jupyter Notebook widget.

    Parameters
    ----------
    target : Hamiltonian or Structure
        The Hamiltonian or crystal Structure object to display.
    initial_rotation : (3, 3) numpy.ndarray or None, optional
        Initial 3x3 view rotation matrix. Default is None.
    initial_distance : float or None, optional
        Initial camera distance from the view origin. Default is None.
    display_options : DisplayOptions or None, optional
        Initial display configuration. Default is None.
    copies : tuple of (int, int, int), optional
        Supercell scaling along each lattice vector. Default is (1, 1, 1).
    rotation_supercell_expansion_max : tuple of (int, int, int), optional
        Maximum denominator along each axis when approximating a
        RotationSupercell propagation vector as a commensurate cell. Default is (10, 10, 10).
    propagation_vector_commensurability_tolerance : float, optional
        Tolerance for matching a propagation vector component to a rational fraction. Default is 1e-4.
    width : int, optional
        Canvas width in pixels. Default is 700.
    height : int, optional
        Canvas height in pixels. Default is 500.

    Returns
    -------
    ViewerJupyter
        An interactive composite widget displaying the 3D crystal structure with controls.

    Raises
    ------
    TypeError
        If target is not an instance of Hamiltonian or Structure.
    """
    if isinstance(target, Structure):
        target = Hamiltonian(target, [])

    if not isinstance(target, Hamiltonian):
        raise TypeError(f"view_jupyter expected Hamiltonian or Structure, got {type(target).__name__}")

    return ViewerJupyter(
        hamiltonian=target,
        initial_rotation=initial_rotation,
        initial_distance=initial_distance,
        display_options=display_options,
        copies=copies,
        rotation_supercell_expansion_max=rotation_supercell_expansion_max,
        propagation_vector_commensurability_tolerance=propagation_vector_commensurability_tolerance,
        width=width,
        height=height,
    )
