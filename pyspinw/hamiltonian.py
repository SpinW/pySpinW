"""Selection of different Hamiltonians"""

from abc import ABC, abstractmethod

import numpy as np

import matplotlib.pyplot as plt
from numpy._typing import ArrayLike

from pyspinw.anisotropy import Anisotropy
from pyspinw.calculations.spinwave import (
    spinwave_calculation as py_spinwave,
    Coupling as PyCoupling,
    MagneticField as PyMagneticField)

from pyspinw.anisotropy import Anisotropy
from pyspinw.cell_offsets import CellOffset
from pyspinw.checks import check_sizes
from pyspinw.coupling import Coupling
from pyspinw.path import Path
from pyspinw.serialisation import SPWSerialisable, SPWSerialisationContext, SPWDeserialisationContext, expects_keys
from pyspinw.structures import Structure
from pyspinw.basis import site_rotations


# pylint: disable=R0903

class Hamiltonian(SPWSerialisable):
    """Hamiltonian base class"""

    serialisation_name = "hamiltonian"

    def __init__(self,
                 structure: Structure,
                 couplings: list[Coupling],
                 anisotropies: list[Anisotropy] | None = None):

        self._structure = structure
        self._couplings = couplings
        self._anisotropies = [] if anisotropies is None else anisotropies

    @property
    def structure(self):
        """ Get the magnetic structure """
        return self._structure

    @property
    def couplings(self):
        """ Get the couplings """
        return self._couplings

    @property
    def anisotropies(self):
        """ Get the anisotropies """
        return self._anisotropies

    @property
    def text_summary(self) -> str:
        """ String giving details of the system """

        lines = ["Sites:"]
        for site in self.structure.sites:
            lines.append(f"  {site}")

        lines.append("Couplings:")
        for exchange in self.couplings:
            lines.append(f"  {exchange}") # vector =", exchange.vector(unit_cell=unit_cell))

        if self.anisotropies:
            lines.append("Anisotropies:")
            for anisotropy in self.anisotropies:
                lines.append(f"  {anisotropy}")

        return "\n".join(lines)

    def print_summary(self):
        """ Print a textual summary to stdout"""
        print(self.text_summary)

    @check_sizes(q_vectors=(-1, 3), field=(3,), allow_nones=True, force_numpy=True)
    def energies_and_intensities(self, q_vectors: np.ndarray, field: ArrayLike | None = None, use_rust: bool=True):
        """Calculate the energy levels of the system for the given q-vectors."""
        # default to Python unless Rust is requested (which it is by default) and available
        coupling_class = PyCoupling
        spinwave_calculation = py_spinwave
        magnetic_field_class = PyMagneticField

        if use_rust:
            try:
                from pyspinw.rust import (
                    spinwave_calculation as rs_spinwave,
                    Coupling as RsCoupling,
                    MagneticField as RsMagneticField)

                coupling_class = RsCoupling
                spinwave_calculation = rs_spinwave
                magnetic_field_class = RsMagneticField


            except ModuleNotFoundError:
                # Silently don't use rust, maybe should give a warning though
                pass

        rust_kw = {'dtype': complex, 'order': 'F'}

        # Get the positions, rotations, moments for the sites
        moments = []
        positions = []
        unique_id_to_index: dict[int, int] = {}
        for index, site in enumerate(self._structure.sites):
            # TODO: Sort out moments for supercells
            moments.append(
                self._structure.supercell.moment(
                    site,
                    cell_offset=CellOffset(0,0,0)
                )
            )

            positions.append(
                self._structure.unit_cell.fractional_to_cartesian(site.ijk))

            unique_id_to_index[site._unique_id] = index

        # Get the field object
        if field is None:
            magnetic_field = None
        else:
            g_tensors = []
            for site in self._structure.sites:
                g_tensors.append(site.g)

            magnetic_field = magnetic_field_class(vector=np.array(field), g_tensors=g_tensors)

        moments = np.array(moments, dtype=float)
        rotations = site_rotations(moments)
        magnitudes = np.sqrt(np.sum(moments**2, axis=1))
        rotations = [rotations[i, :, :] for i in range(rotations.shape[0])]

        # Convert the couplings
        couplings: list[Coupling] = []
        for input_coupling in self._couplings:
            # Normal coupling

            coupling = coupling_class(
                unique_id_to_index[input_coupling.site_1._unique_id],
                unique_id_to_index[input_coupling.site_2._unique_id],
                np.array(input_coupling.coupling_matrix, **rust_kw),
                input_coupling.vector(self.structure.unit_cell)
            )

            couplings.append(coupling)

            # Reversed coupling

            coupling = coupling_class(
                unique_id_to_index[input_coupling.site_2._unique_id],
                unique_id_to_index[input_coupling.site_1._unique_id],
                np.array(input_coupling.coupling_matrix.T, **rust_kw),
                -input_coupling.vector(self.structure.unit_cell)
            )

            couplings.append(coupling)

        # Add in anisotropies as spinwave_calculation couplings
        for input_anisotropy in self._anisotropies:

            anisotropy = coupling_class(
                unique_id_to_index[input_anisotropy.site._unique_id],
                unique_id_to_index[input_anisotropy.site._unique_id],
                np.array(input_anisotropy.anisotropy_matrix.T, **rust_kw),
                inter_site_vector=np.array([0,0,0], dtype=float)
            )

            couplings.append(anisotropy)

        result = spinwave_calculation(
                        rotations=rotations,
                        magnitudes=magnitudes,
                        q_vectors=q_vectors,
                        couplings=couplings,
                        positions=positions,
                        field=magnetic_field)

        return result[0], result[1]

    def sorted_positive_energies(self,
                                 path: Path,
                                 field: ArrayLike | None = None,
                                 use_rust: bool = True) -> list[np.ndarray]:

        """ Return energies as series corresponding to q, sorted by energy """
        energy, _ = self.energies_and_intensities(path.q_points(), field=field, use_rust=use_rust)

        # Sort the energies
        energy = np.sort(energy.real, axis=1)

        # return the top half (positive)
        n_energies = energy.shape[1]
        energies = [energy[:, n_energies - i - 1] for i in range(n_energies//2)]

        return energies

    def energy_plot(self,
                    path: Path,
                    field: ArrayLike | None = None,
                    show: bool=True,
                    new_figure: bool=True,
                    use_rust: bool=True):

        """ Create a spaghetti diagram """

        if new_figure:
            plt.figure("Energy")


        x_values = path.x_values()

        for series in self.sorted_positive_energies(path, field=field, use_rust=use_rust):
            plt.plot(x_values, series, 'k')

        path.format_plot(plt)

        if show:
            plt.show()

    def _serialise(self, context: SPWSerialisationContext) -> dict:
        return {"magnetic_structure": self.structure._serialise(context),
                "couplings": [coupling._serialise(context) for coupling in self.couplings],
                "anisotropies": [anisotropy._serialise(context) for anisotropy in self.anisotopies]}

    @staticmethod
    @expects_keys("magnetic_structure, couplings, anisotropies")
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        structure = Structure._deserialise(json["magnetic_structure"], context)
        couplings = [Coupling._deserialise(coupling, context) for coupling in json["couplings"]]
        anisotropies = [Anisotropy._deserialise(anisotropy, context) for anisotropy in json["anisotropies"]]

        return Hamiltonian(structure, couplings, anisotropies)