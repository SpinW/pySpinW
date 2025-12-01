"""Selection of different Hamiltonians"""

from abc import ABC, abstractmethod

import numpy as np

import matplotlib.pyplot as plt

from pyspinw.calculations.spinwave import spinwave_calculation as py_spinwave, Coupling as PyCoupling

from pyspinw.cell_offsets import CellOffset
from pyspinw.coupling import Coupling
from pyspinw.path import Path
from pyspinw.serialisation import SPWSerialisable
from pyspinw.structures import Structure
from pyspinw.basis import find_aligned_basis, site_rotations


# pylint: disable=R0903

class Hamiltonian(SPWSerialisable):
    """Hamiltonian base class"""

    serialisation_name = "hamiltonian"

    def __init__(self,
                 structure: Structure,
                 couplings: list[Coupling]):

        self._structure = structure
        self._couplings = couplings

    @property
    def structure(self):
        """ Get the magnetic structure """
        return self._structure

    @property
    def couplings(self):
        """ Get the couplings """
        return self._couplings

    def energies(self, q_vectors: np.ndarray, use_rust: bool=True):
        """Calculate the energy levels of the system for the given q-vectors."""
        # default to Python unless Rust is requested (which it is by default) and available
        coupling_class = PyCoupling
        spinwave_calculation = py_spinwave

        if use_rust:
            try:
                from pyspinw.rust import spinwave_calculation as rs_spinwave, Coupling as RsCoupling

                coupling_class = RsCoupling
                spinwave_calculation = rs_spinwave


            except ModuleNotFoundError:
                # Silently don't use rust, maybe should give a warning though
                pass

        rust_kw = {'dtype': complex, 'order': 'F'}

        # Get the rotations for the sites
        moments = []
        unique_id_to_index: dict[int, int] = {}
        for index, site in enumerate(self._structure.sites):
            # TODO: Sort out moments for supercells
            moments.append(
                self._structure.supercell.moment(
                    site,
                    cell_offset=CellOffset(0,0,0)
                )
            )

            unique_id_to_index[site._unique_id] = index

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

        energies = spinwave_calculation(rotations, magnitudes, q_vectors, couplings)

        return energies


    def energy_plot(self, path: Path, show=True, new_figure=True):
        """ Create a spaghetti diagram """
        energy = self.energies(path.q_points())

        # Sort the energies
        energy = np.sort(energy.real, axis=1)

        if new_figure:
            plt.figure("Energy")

        n_energies = energy.shape[1]

        x_values = path.x_values()
        for i in range(n_energies//2):
            plt.plot(x_values, energy[:, n_energies-i-1], 'k')

        path.format_plot(plt)

        if show:
            plt.show()

    def _serialise(self) -> dict:
        return {"magnetic_structure": self.structure._serialise(),
                "couplings": [coupling._serialise() for coupling in self.couplings]}

