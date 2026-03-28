"""Selection of different Hamiltonians"""
import logging
from abc import ABC, abstractmethod
import re
from collections import Counter

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
from pyspinw.site import LatticeSite
from pyspinw.structures import Structure
from pyspinw.basis import site_rotations
from pyspinw.symmetry.supercell import TrivialSupercell, RotationSupercell

# pylint: disable=R0903

logger = logging.Logger("pyspinw.hamiltonian")

def omegasum(energy: ArrayLike, intensity: ArrayLike, tol: float=1e-5, zeroint: int=0):
    """ Removes degenerate and ghost (zero-intensity) modes from spectrum """
    energy, intensity = (np.array(energy), np.array(intensity))
    en_out, int_out = (energy * np.nan, intensity * np.nan)
    if tol > 0:
        energy = np.round(energy / tol) * tol
    if zeroint > 0:
        energy[np.where(np.abs(np.real(intensity)) < zeroint)] = np.nan
    for iQ in range(energy.shape[0]):
        eu = np.unique(energy[iQ,:])
        en_out[iQ,:len(eu)] = eu
        for iE in range(len(eu)):
            int_out[iQ, iE] = np.sum(intensity[iQ, np.where(np.abs(energy[iQ,:] - en_out[iQ, iE]) < tol)])
    return en_out, int_out


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

    def couplings_by_name(self, regex):
        """ Get list of couplings whose names match the regex """
        return [coupling for coupling in self.couplings if re.match(regex, coupling.name) is not None]


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

    def _expand_with_mapping(self) -> tuple["Hamiltonian",
                                            dict[tuple[int, tuple[int, int, int]], LatticeSite],
                                            list[int],
                                            list[int]]:
        """ Expand the supercell structure into a single cell structure and return the mapping between hamiltonians

        This should only be used internally, and its not very user friendly

        The site mapping is a dict from (original site UID, offset) -> LatticeSite
        The coupling mapping is a len(new couplings) list of indices for the original couplings
        The anisotropy mapping is a len(new anisotropies) list of indices for the original anisotropies
        """
        bigger_cell, site_mapping = self.structure._expansion_site_mapping()

        new_couplings = []
        new_anisotropies = []

        si, sj, sk = self.structure.supercell.cell_size()

        coupling_mapping = []
        anisotropy_mapping = []

        for first_site_offset in self.structure.supercell.cells():
            for original_index, coupling in enumerate(self.couplings):
                # Convert the offset in the coupling, into
                #  1) an offset in the supercell, and
                #  2) an offset between supercells
                # basically just a divmod

                # Convert offsets into "absolute offsets"

                second_site_offset = coupling.cell_offset.vector + first_site_offset.vector

                oi, oj, ok = second_site_offset

                ni, ri = divmod(oi, si)
                nj, rj = divmod(oj, sj)
                nk, rk = divmod(ok, sk)

                new_cell_offset = CellOffset(ni, nj, nk)
                second_site_lookup_value = (ri, rj, rk)

                target_site_1 = site_mapping[(coupling.site_1.unique_id, first_site_offset.as_tuple)]
                target_site_2 = site_mapping[(coupling.site_2.unique_id, second_site_lookup_value)]

                # Create the new coupling using their update method, which copies everything not specified
                new_couplings.append(
                    coupling.updated(
                        site_1=target_site_1,
                        site_2=target_site_2,
                        cell_offset=new_cell_offset))

                # Add entry to mapping
                coupling_mapping.append(original_index)

            for original_index, anisotropy in enumerate(self.anisotropies):
                target_site = site_mapping[(anisotropy.site.unique_id, first_site_offset.as_tuple)]

                # Update-copy anisotropy
                new_anisotropies.append(anisotropy.updated(site=target_site))

                # Add entry to mapping
                anisotropy_mapping.append(original_index)

        structure = Structure(
            sites=[site for site in site_mapping.values()],
            unit_cell=bigger_cell,
            spacegroup=self.structure.spacegroup.for_supercell(self.structure.supercell),
            supercell=TrivialSupercell(scaling=(1,1,1))
        )

        return (Hamiltonian(structure=structure, couplings=new_couplings, anisotropies=new_anisotropies),
                site_mapping, coupling_mapping, anisotropy_mapping)

    def expanded(self):
        """ Expand the supercell structure into a single cell structure """
        expanded, _, _, _ = self._expand_with_mapping()
        return expanded


    def print_summary(self):
        """ Print a textual summary to stdout"""
        print(self.text_summary)

    @check_sizes(q_vectors=(-1, 3), field=(3,), allow_nones=True, force_numpy=True)
    def energies_and_intensities(self,
                                 q_vectors: np.ndarray,
                                 field: ArrayLike | None = None,
                                 use_rust: bool=True,
                                 use_rotating: bool=True,
                                 ):
        """Calculate the energy levels of the system for the given q-vectors."""
        #
        # Set up choice of calculation
        #

        # Rotating frame calculations should only be run on RotationSupercells
        use_rotating = use_rotating and isinstance(self.structure.supercell, RotationSupercell)

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
                logger.warning("Failed to load rust core, falling back to python")

        else:
            logger.info("Using Python core code")

        rust_kw = {'dtype': complex, 'order': 'F'}

        #
        # Set up the system
        #
        if all([hasattr(self.structure.supercell, v) for v in ['propagation_vector', 'perpendicular']]):
            if use_rotating:
                expanded, scaling = (self, 1.)
                nvec = self.structure.supercell.perpendicular
                rotating_frame = [self.structure.supercell.propagation_vector._vector, nvec / np.linalg.norm(nvec)]
            else:
                newstruc = Structure(**{k:getattr(self.structure, k) for k in ['sites', 'unit_cell', 'spacegroup']},
                               supercell=self.structure.supercell.approximant())
                expanded = Hamiltonian(newstruc, self.couplings, self.anisotropies).expanded()
                scaling, rotating_frame = (newstruc.supercell.scaling, None)
        else:
            if use_rotating:
                logger.warning("Cannot do rotating frame calculation propagation vector or plane normal not specified")
            expanded, scaling, rotating_frame = (self.expanded(), self.structure.supercell.scaling, None)

        # Get the positions, rotations, moments for the sites
        moments = []
        positions = []
        unique_id_to_index: dict[int, int] = {}
        for index, site in enumerate(expanded.structure.sites):
            # TODO: Sort out moments for supercells
            moments.append(site.base_moment)

            positions.append(site.ijk)

            unique_id_to_index[site._unique_id] = index

        # Get the field object
        if field is None:
            magnetic_field = None
        else:
            g_tensors = []
            for site in expanded.structure.sites:
                g_tensors.append(site.g)
            magnetic_field = magnetic_field_class(
                                vector=np.array(field, **rust_kw),
                                g_tensors=np.array(g_tensors, **rust_kw))

        moments = np.array(moments, dtype=float)
        rotations = site_rotations(moments)
        magnitudes = np.sqrt(np.sum(moments**2, axis=1))
        rotations = np.array([rotations[i, :, :] for i in range(rotations.shape[0])], **rust_kw)

        # Convert the couplings
        couplings: list[Coupling] = []
        for input_coupling in expanded.couplings:
            # Normal coupling

            coupling = coupling_class(
                unique_id_to_index[input_coupling.site_1._unique_id],
                unique_id_to_index[input_coupling.site_2._unique_id],
                np.array(input_coupling.coupling_matrix, **rust_kw),
                input_coupling.cell_offset.vector.astype('double')
            )

            couplings.append(coupling)

            # Reversed coupling

            coupling = coupling_class(
                unique_id_to_index[input_coupling.site_2._unique_id],
                unique_id_to_index[input_coupling.site_1._unique_id],
                np.array(input_coupling.coupling_matrix.T, **rust_kw),
                -input_coupling.cell_offset.vector.astype('double')
            )

            couplings.append(coupling)

        ## This shouldn't be needed now things are (probably) fixed - keeping in case
        # Remove duplicate couplings
        # couplings = [c1 for ic, c1 in enumerate(couplings)
        #              if all([c1 != c2 for c2 in couplings[ic+1:]])]

        # Add in anisotropies as spinwave_calculation couplings
        for input_anisotropy in expanded.anisotropies:

            anisotropy = coupling_class(
                unique_id_to_index[input_anisotropy.site._unique_id],
                unique_id_to_index[input_anisotropy.site._unique_id],
                # Factor 2 here is to agree with Matlab code (need to check if is correct)
                np.array(input_anisotropy.anisotropy_matrix.T, **rust_kw) * 2.,
                inter_site_vector=np.array([0,0,0], dtype=float)
            )

            couplings.append(anisotropy)

        result = spinwave_calculation(
                        rotations=rotations,
                        magnitudes=magnitudes,
                        q_vectors=q_vectors * scaling,
                        couplings=couplings,
                        positions=positions,
                        rlu_to_cart=np.linalg.inv(self.structure.unit_cell._xyz).T * 2 * np.pi,
                        field=magnetic_field,
                        rotating_frame=rotating_frame)

        # Applies a rescaling to agree with Matlab code for Sab
        # Toth & Lake eq (46) gives a 1/(2Natom) prefactor but the Matlab code uses 1/(2*Ncell)
        scale_factor = rotations.shape[0] / np.prod(scaling)
        intensity = [res * scale_factor for res in result[1]]

        return result[0], intensity

    def spaghetti_plot(self,
             path: Path,
             field: ArrayLike | None = None,
             show: bool=True,
             new_figure: bool=True,
             use_rust: bool=True,
             use_rotating: bool=False,
             scale: str='linear'):
        """ Create a spaghetti diagram with energy top and intensity bottom """
        if new_figure:
            fig, axs = plt.subplots(2, 1)
        else:
            fig = plt.gcf()
            axs = fig.get_axes()
            for ii in range(len(axs), 2):
                axs.append(fig.add_subplot(2,1,ii+1))

        x_values = path.x_values()
        energy, intensity = omegasum(*self.energies_and_intensities(path.q_points(),
                                     field=field, use_rust=use_rust, use_rotating=use_rotating))
        n_mode = energy.shape[1]
        for series in zip(*([v[:, n_mode - i - 1] for i in range(n_mode)] for v in (energy, intensity))):
            axs[0].plot(x_values, series[0], 'k')
            axs[1].plot(x_values, series[1])
        if 'log' in scale:
            axs[1].set_yscale('log')

        path.format_plot(axs[0])
        path.format_plot(axs[1])

        if show:
            plt.show()
        else:
            return fig

    def parameterize(self, *parameters: tuple[Coupling, str] | tuple[str, str] | str) \
            -> "HamiltonianParameterization":
        """ Get a function that maps floats to a hamiltonian with the floats controlling the specified parameters"""

        return HamiltonianParameterization(self, *parameters)

    def sorted_positive_energies(self,
                                 path: Path,
                                 field: ArrayLike | None = None,
                                 use_rust: bool = True) -> list[np.ndarray]:
        """ Return energies as series corresponding to q, sorted by energy """
        energy, intensities = self.energies_and_intensities(path.q_points(), field=field, use_rust=use_rust)

        energy = np.array(energy)

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


class HamiltonianParameterization:
    def __init__(self,
                 hamiltonian: Hamiltonian,
                 *parameters: list[tuple[Coupling, str] | tuple[str, str] | str]):

        self.hamiltonian = hamiltonian

        # Create a list of parameters that will be updated
        base_parameter_definitions = []
        for parameter_data in parameters:

            # input case: str, split into a tuple, then parse like others
            if isinstance(parameter_data, str):

                parts = parameter_data.split(".")
                if len(parts) != 2:
                    raise ValueError("Expected parameter definition to be of the form 'coupling_name.parameter', "
                                     f"got '{parameter_data}'")

                parameter_data = (parts[0], parts[1])

            # Deal with tuples
            if isinstance(parameter_data, tuple):

                if len(parameter_data) != 2:
                    raise ValueError(f"Expected tuple entries to be length 2, got {parameter_data}")

                coupling, parameter = parameter_data
                if isinstance(coupling, Coupling):
                    # input case: (Coupling, str)

                    base_parameter_definitions.append([(coupling, parameter)])

                elif isinstance(coupling, str):
                    # input case: (str, str)

                    couplings = hamiltonian.couplings_by_name(coupling)

                    if len(couplings) == 0:
                        logger.warning(f"Coupling regex {coupling} does not match any couplings")

                    elif len(couplings) > 1:
                        logger.warning(f"Coupling regex {coupling} matches multiple couplings, adding all")

                    base_parameter_definitions.append([(coupling, parameter) for coupling in couplings])

                else:
                    raise TypeError(f"Expected first component of tuple to be Coupling or str, got {type(coupling)}")

            else:
                # input case: bad
                raise TypeError(f"Expected parameters to be tuples of Exchange/str and str, "
                                f"or str, got {type(parameter_data)}")

        # Check that there is no conflicts, this means that each parameter is only set by only one entry
        # Basically, we can just check for duplicates

        multicounts = [item for item, count in Counter(sum(base_parameter_definitions, [])).items() if count > 1]

        if len(multicounts) > 0:
            raise ValueError(f"Multiple parameters assigned to the same parameter {multicounts}")

        # Check that the couplings have the required parameters
        for parameter_definition in base_parameter_definitions:
            for coupling, attribute in parameter_definition:
                if attribute not in coupling.parameters:
                    valid_parameters = ", ".join(coupling.parameters)
                    raise TypeError(f"{coupling} does not have parameter '{attribute}', it has: {valid_parameters}")


        # Convert the coupling to index in the list of couplings
        unique_id_to_index = {coupling.unique_id: index for index, coupling in enumerate(hamiltonian.couplings)}
        self.parameter_definitions: list[list[tuple[int, str]]] = []

        for parameter_definition in base_parameter_definitions:
            indexed_parameter_definition: list[tuple[int, str]] = []
            for coupling, attribute in parameter_definition:
                index = unique_id_to_index[coupling.unique_id]
                indexed_parameter_definition.append((index, attribute))
            self.parameter_definitions.append(indexed_parameter_definition)

        # Number is useful
        self.n_parameters = len(self.parameter_definitions)

    def energy_plot(self,
                    parameter_values: ArrayLike,
                    path: Path,
                    field: ArrayLike | None = None,
                    show: bool = True,
                    colormap_name: str = 'jet',
                    new_figure: bool = True,
                    use_rust: bool = True):
        """ Show a plot of the energies """

        # Check/regularise input values
        parameter_values = np.array(parameter_values)
        if len(parameter_values.shape) == 1:
            parameter_values = parameter_values.reshape(-1, 1)

        if parameter_values.shape[1] != self.n_parameters:
            raise ValueError("Second dimension size of parameter_values should match number of parameters")

        n_curves = parameter_values.shape[0]

        #
        # Do the plotting
        #

        # Colours
        try:
            cmap = plt.get_cmap(colormap_name)
        except ValueError as ve:
            logger.warning(str(ve))
            cmap = plt.get_cmap('jet')

        colors = cmap(np.linspace(0, 1, n_curves+2)[1:-1])

        # Figure
        if new_figure:
            plt.figure("Energy")

        x_values = path.x_values()

        for i in range(n_curves):
            ham = self(*parameter_values[i, :])
            for series in ham.sorted_positive_energies(path, field=field, use_rust=use_rust):
                plt.plot(x_values, series, color=colors[i])

        path.format_plot(plt)

        if show:
            plt.show()



    def __call__(self, *parameters: float) -> Hamiltonian:
        """ Get the Hamiltonian with parameters set """
        if len(parameters) != self.n_parameters:
            raise ValueError(f"Expected {self.n_parameters} parameters, got {len(parameters)}")

        new_couplings = [coupling for coupling in self.hamiltonian.couplings]
        for parameter_definition, value in zip(self.parameter_definitions, parameters):
            for (coupling_index, attribute) in parameter_definition:
                new_couplings[coupling_index] = new_couplings[coupling_index].updated(**{attribute: value})


        return Hamiltonian(self.hamiltonian.structure, new_couplings, self.hamiltonian.anisotropies)