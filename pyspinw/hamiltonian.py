"""Selection of different Hamiltonians"""
import logging
from collections.abc import Callable
import re
from collections import Counter
from typing import Sequence, Union

import numpy as np

import matplotlib.pyplot as plt
from numpy._typing import ArrayLike

from pyspinw.calculations.energy_minimisation import ClassicalEnergyMinimisation, Free, \
    InitialRandomisation, Fixed, Planar
from pyspinw.calculations.spinwave import (
    spinwave_calculation as py_spinwave,
    Coupling as PyCoupling,
    MagneticField as PyMagneticField)

from pyspinw.anisotropy import Anisotropy
from pyspinw.cell_offsets import CellOffset
from pyspinw.checks import check_sizes
from pyspinw.exchange import Exchange
from pyspinw.path import Path
from pyspinw.serialisation import SPWSerialisable, SPWSerialisationContext, SPWDeserialisationContext, expects_keys
from pyspinw.site import LatticeSite
from pyspinw.structures import Structure
from pyspinw.basis import site_rotations
from pyspinw.symmetry.supercell import TiledSupercell, RotationSupercell
from pyspinw.units import IntensityUnits, intensity_units

# pylint: disable=R0903

logger = logging.Logger("pyspinw.hamiltonian")

def uniquetol(values: ArrayLike, tol: float=1e-5):
    """ Returns floating point unique values within a given tolerance """
    v = np.sort(values)
    idif = np.append(True, np.diff(v))
    return v[idif > tol]

def omegasum(energy: ArrayLike, intensity: ArrayLike, tol: float=1e-5, zeroint: int=0, is_series: bool=True):
    """ Removes degenerate and ghost (zero-intensity) modes from spectrum """
    energy, intensity = (np.array(energy), np.array(intensity))
    en_out, int_out = (energy * np.nan, intensity * np.nan)
    if tol > 0:
        energy = np.round(energy / tol) * tol
    if zeroint > 0:
        energy[np.where(np.abs(np.real(intensity)) < zeroint)] = np.nan
    for iQ in range(energy.shape[0]):
        eu = uniquetol(energy[iQ,:])
        en_out[iQ,:len(eu)] = eu
        for iE in range(len(eu)):
            int_out[iQ, iE] = np.sum(intensity[iQ, np.where(np.abs(energy[iQ,:] - en_out[iQ, iE]) < tol)])
    nans = np.isnan(en_out)
    nanmodes = np.sum(nans, axis=0)
    if is_series:
        # Ensure there are the same number of modes throughout
        for col in set([int(idx) for accidentals in np.where((nanmodes > 0) * (nanmodes < en_out.shape[0]))[0]
                        for idx in np.where(nans[:,accidentals])[0]]):
            idnxt = 1 if col < en_out.shape[0]-1 else -1
            idnan = np.where(~np.isnan(en_out[col + idnxt,:]))[0]
            idx = np.array([np.nanargmin(np.abs(en_out[col + idnxt,i] - en_out[col])) for i in idnan])
            en_out[col,idnan] = np.take(en_out[col, :], idx)
            int_out[col,idnan] = np.take(int_out[col, np.where(~np.isnan(int_out[col,:]))] / np.bincount(idx), idx)
    # Removes columns which are all NaNs
    en_out = np.delete(en_out, np.where(nanmodes == en_out.shape[0])[0], axis=1)
    int_out = np.delete(int_out, np.where(nanmodes == en_out.shape[0])[0], axis=1)
    return en_out, int_out

def egrid(energy: ArrayLike, intensity: ArrayLike, evect: ArrayLike, dE: float | Callable):
    """Bins a set of energy/intensity into a spectrum with Gaussian broadening"""
    energy, intensity = (np.array(energy), np.array(intensity))
    esigma = dE(evect) if isinstance(dE, Callable) else np.ones(evect.shape) * dE
    esigma /= 2 * np.sqrt(2* np.log(2))  # Convert from FWHM
    spec = np.zeros((energy.shape[0], len(evect)))
    ew_norm = np.diff(evect, append=evect[-1] + (evect[-1] - evect[-2])) / (np.sqrt(2 * np.pi) * esigma)
    for iE in range(energy.shape[1]):
        spec += intensity[:,iE,np.newaxis] * ew_norm * \
            np.exp(-0.5 * ((evect[np.newaxis,:] - energy[:,iE,np.newaxis]) / esigma)**2)
    spec[np.where(spec < np.finfo(np.float32).eps)] = np.nan
    return spec

ParametrizationType = Union[str,
                       list[str],
                       list[tuple[Union[Exchange, Anisotropy, str], str]],
                       tuple[Sequence[Union[Exchange, Anisotropy, str]], str],
                       tuple[Exchange | Anisotropy | str, str]]

def _regularise_parameters(hamiltonian: "Hamiltonian", parameter_data: ParametrizationType):
    """ Convert many options for defining parameters into a more regular form

    This is used in the parameterisation method/class
    """
    # Step 1: convert everything to a Sequence

    if isinstance(parameter_data, str):
        parameter_data = [parameter_data]

    elif isinstance(parameter_data, tuple):
        if len(parameter_data) != 2:
            raise ValueError(f"Expected tuple parameter definition to have length 2 (got {parameter_data})")

        if not isinstance(parameter_data[1], str):
            raise ValueError(f"Expected second component of tuple to be (got {type(parameter_data)}, {parameter_data})")

        if isinstance(parameter_data[0], Sequence):
            parameter_data = [(datum, parameter_data[1]) for datum in parameter_data[0]]

        elif isinstance(parameter_data[0], (Exchange, Anisotropy, str)):
            parameter_data = [parameter_data]

        else:
            raise TypeError(f"Expected parameters to be defined by Exchange, Anisotropy or str,"
                            f" got {type(parameter_data[0])}: {parameter_data[0]}")

    elif isinstance(parameter_data, list):
        pass

    else:
        raise TypeError(f"Expected to parameter definitions type to be str, Sequence, or tuple, "
                        f"got {type(parameter_data)}: {parameter_data}")

    #
    # Step 2: convert plain strings into tuples
    #

    new_parameter_data = []
    for datum in parameter_data:
        if isinstance(datum, str):
            parts = datum.split(".")

            if len(parts) != 2:
                raise ValueError(f"Expected strings to be of the form exchange_name.parameter_name, got '{datum}'")

            new_parameter_data.append((parts[0], parts[1]))

        else:
            new_parameter_data.append(datum)

    parameter_data = new_parameter_data

    #
    # Step 3: Check second part is always a string, and that tuples are the right length
    #

    for datum in parameter_data:
        if not isinstance(datum, tuple):
            raise Exception("A previous error has not been caught, parameters should all be tuples")

        if len(datum) != 2:
            raise ValueError(f"Expected all tuple entries to be length 2, found {datum}")

        if not isinstance(datum[1], str):
            raise ValueError(f"Expected second part of tuple to be a string, got {type(datum[1])}: {datum[1]}")

    #
    # Step 4: Resolve any string names, everything should now be tuple2
    #

    new_parameter_data = []
    for target, parameter in parameter_data:
        if isinstance(target, str):

            exchanges = hamiltonian.exchanges_by_name(target)

            if len(exchanges) == 0:
                raise ValueError(f"Could not find exchanges that match {target}")

            if len(exchanges) > 1:
                logger.warning(f"Found multiple matches for '{target}', attempting to use them all: {exchanges}")

            for exchange in exchanges:
                new_parameter_data.append((exchange, parameter))

        else:
            new_parameter_data.append((target, parameter))

    parameter_data = new_parameter_data

    #
    # Step 5: Split into exchanges and anisotopies
    #

    exchanges = []
    anisotropies = []
    for target, parameter in parameter_data:
        if isinstance(target, Exchange):
            exchanges.append((target, parameter))
        elif isinstance(target, Anisotropy):
            anisotropies.append((target, parameter))

    return exchanges, anisotropies


class Hamiltonian(SPWSerialisable):
    """Hamiltonian base class"""

    serialisation_name = "hamiltonian"

    def __init__(self,
                 structure: Structure,
                 exchanges: list[Exchange],
                 anisotropies: list[Anisotropy] | None = None):

        self._structure = structure
        self._exchanges = exchanges
        self._anisotropies = [] if anisotropies is None else anisotropies

    @property
    def structure(self):
        """ Get the magnetic structure """
        return self._structure

    @property
    def exchanges(self):
        """ Get the exchanges """
        return self._exchanges

    def exchanges_by_name(self, regex):
        """ Get list of exchanges whose names match the regex """
        return [exchange for exchange in self.exchanges if re.match(regex, exchange.name) is not None]


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

        lines.append("Exchanges:")
        for exchange in self.exchanges:
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
        The exchange mapping is a len(new exchanges) list of indices for the original exchanges
        The anisotropy mapping is a len(new anisotropies) list of indices for the original anisotropies
        """
        bigger_cell, site_mapping = self.structure._expansion_site_mapping()

        new_exchanges = []
        new_anisotropies = []

        si, sj, sk = self.structure.supercell.cell_size()

        exchange_mapping = []
        anisotropy_mapping = []

        for first_site_offset in self.structure.supercell.cells():
            for original_index, exchange in enumerate(self.exchanges):
                # Convert the offset in the exchange, into
                #  1) an offset in the supercell, and
                #  2) an offset between supercells
                # basically just a divmod

                # Convert offsets into "absolute offsets"

                second_site_offset = exchange.cell_offset.vector + first_site_offset.vector

                oi, oj, ok = second_site_offset

                ni, ri = divmod(oi, si)
                nj, rj = divmod(oj, sj)
                nk, rk = divmod(ok, sk)

                new_cell_offset = CellOffset(ni, nj, nk)
                second_site_lookup_value = (ri, rj, rk)

                target_site_1 = site_mapping[(exchange.site_1.unique_id, first_site_offset.as_tuple)]
                target_site_2 = site_mapping[(exchange.site_2.unique_id, second_site_lookup_value)]

                # Create the new exchange using their update method, which copies everything not specified
                new_exchanges.append(
                    exchange.updated(
                        site_1=target_site_1,
                        site_2=target_site_2,
                        cell_offset=new_cell_offset))

                # Add entry to mapping
                exchange_mapping.append(original_index)

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
            supercell=TiledSupercell(scaling=(1, 1, 1))
        )

        return (Hamiltonian(structure=structure, exchanges=new_exchanges, anisotropies=new_anisotropies),
                site_mapping, exchange_mapping, anisotropy_mapping)

    def expanded(self):
        """ Expand the supercell structure into a single cell structure """
        expanded, _, _, _ = self._expand_with_mapping()
        return expanded


    def sites_by_name(self, regex) -> list[LatticeSite]:
        """ Get sites where name matches regex"""
        return self.structure.sites_by_name(regex)

    def print_summary(self):
        """ Print a textual summary to stdout"""
        print(self.text_summary)

    @check_sizes(q_vectors=(-1, 3), field=(3,), allow_nones=True, force_numpy=True)
    def energies_and_intensities(self,
                                 q_vectors: np.ndarray,
                                 field: ArrayLike | None = None,
                                 use_rust: bool=True,
                                 use_rotating: bool=True,
                                 intensity_unit: IntensityUnits | str='cell',
                                 ):
        """Calculate the energy levels of the system for the given q-vectors.

        :param q_vectors: *required* An array of q-vectors
        :param field: Optional field direction
        :param use_rust: Whether to use Rust or Python calculator (default: True)
        :param use_rotating: Whether to use the rotating frame calculator if possible (default: True)
        :param intensity_unit: Whether to normalise intensity per unit cell or spin (default: 'cell')
        """
        intensity_unit = intensity_units(intensity_unit)
        if intensity_unit == IntensityUnits.BARNPERATOM or intensity_unit == IntensityUnits.BARNPERCELL:
            raise NotImplementedError('Intensity scale in barn/sr/meV/atom not yet implemented.')

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
                expanded = Hamiltonian(newstruc, self.exchanges, self.anisotropies).expanded()
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
            moments.append(site.base_spin)

            positions.append(site.ijk)

            unique_id_to_index[site._unique_id] = index

        # Get the field object
        if field is None:
            magnetic_field = None
        else:
            g_tensors = []
            for site in expanded.structure.sites:
                g_tensors.append(site.g)
            # Factor 2 in field below to agree with Matlab code (like with SIA term in L474)
            magnetic_field = magnetic_field_class(
                                vector=2.0 * np.array(field, **rust_kw),
                                g_tensors=np.array(g_tensors, **rust_kw))

        moments = np.array(moments, dtype=float)
        rotations = site_rotations(moments)
        magnitudes = np.sqrt(np.sum(moments**2, axis=1))
        rotations = np.array([rotations[i, :, :] for i in range(rotations.shape[0])], **rust_kw)

        # Convert the couplings
        couplings: list[Exchange] = []
        for input_exchange in expanded.exchanges:
            # Normal coupling

            coupling = coupling_class(
                unique_id_to_index[input_exchange.site_1._unique_id],
                unique_id_to_index[input_exchange.site_2._unique_id],
                np.array(input_exchange.exchange_matrix, **rust_kw),
                input_exchange.cell_offset.vector.astype('double')
            )

            couplings.append(coupling)

            # Reversed coupling

            coupling = coupling_class(
                unique_id_to_index[input_exchange.site_2._unique_id],
                unique_id_to_index[input_exchange.site_1._unique_id],
                np.array(input_exchange.exchange_matrix.T, **rust_kw),
                -input_exchange.cell_offset.vector.astype('double')
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
                        rlu_to_cart=np.linalg.inv(expanded.structure.unit_cell._xyz).T * 2 * np.pi,
                        field=magnetic_field,
                        rotating_frame=rotating_frame)

        # Applies a rescaling to agree with Matlab code for Sab
        # Toth & Lake eq (46) gives a 1/(2Natom) prefactor but the Matlab code uses 1/(2*Ncell)
        if intensity_unit == IntensityUnits.PERCELL:
            scale_factor = rotations.shape[0] / np.prod(scaling)
            intensity = [res * scale_factor for res in result[1]]
        else:
            intensity = result[1]

        return result[0], intensity

    def spaghetti_plot_dual(self,
                            path: Path,
                            field: ArrayLike | None = None,
                            show: bool=True,
                            new_figure: bool=True,
                            use_rust: bool=True,
                            use_rotating: bool=True,
                            intensity_unit: IntensityUnits | str = 'cell',
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
        energy, intensity = omegasum(*self.energies_and_intensities(path.q_points(), field,
                                                                    use_rust, use_rotating, intensity_unit))
        n_mode = energy.shape[1]
        for series in zip(*([v[:, n_mode - i - 1] for i in range(n_mode)] for v in (energy, intensity))):
            axs[0].plot(x_values, series[0], 'k')
            axs[1].plot(x_values, series[1], 'k')
        if 'log' in scale:
            axs[1].set_yscale('log')
        axs[0].set_ylabel('Magnon Energy (meV)')
        axs[1].set_ylabel('Magnon Intensity')

        path.format_plot(axs[0])
        path.format_plot(axs[1])

        if show:
            plt.show()
        else:
            return fig

    def spaghetti_plot(self,
                       path: Path,
                       evect: ArrayLike | None = None,
                       dE: float | Callable | None = None,
                       vmin: float=0,
                       vmax: float | None = None,
                       field: ArrayLike | None = None,
                       show: bool=True,
                       new_figure: bool=True,
                       use_rust: bool=True,
                       use_rotating: bool=True,
                       intensity_unit: IntensityUnits | str = 'cell',
                       scale: str='linear'):
        """ Create a spaghetti diagram with intensity as colorfill overplotted by mode energies """
        if new_figure:
            fig, ax = plt.subplots()
        else:
            fig = plt.gcf()
            ax = fig.get_axes()[0]

        x_values = path.x_values()
        energy, intensity = omegasum(*self.energies_and_intensities(path.q_points(), field,
                                                                    use_rust, use_rotating, intensity_unit))
        if evect is None:
            emax = np.nanmax(energy)
            denom = 10**np.floor(np.log10(emax))
            evect = np.linspace(0, np.ceil(emax / denom) * denom, 100)
        if dE is None:
            dE = np.mean(np.diff(evect)) * 2.35
        spec = egrid(energy, intensity, evect, dE)
        if vmax is None:
            vmax = np.nanmax(spec[:, np.where(evect > max(np.max(evect)/10, 0.1))]) / 10.
        mesh = ax.pcolormesh(x_values, evect, spec.T, vmin=vmin, vmax=vmax)
        ax.plot(x_values, energy, 'k')
        ax.set_ylim(0, np.max(evect))
        fig.colorbar(mesh, ax=ax)
        path.format_plot(ax)
        ax.set_ylabel('Magnon Energy (meV)')

        if show:
            plt.show()
        else:
            return fig

    def parameterize(self,
                     *parameters: ParametrizationType,
                     find_ground_state_with: dict | None = None) -> "HamiltonianParameterization":
        """ Get a function that maps floats to a hamiltonian with the floats controlling the specified parameters"""
        return HamiltonianParameterization(self,*parameters,
                                           find_ground_state_with=find_ground_state_with)

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

    def ground_state(self,
                      fixed: list[LatticeSite] | None = None,
                      planar: list[LatticeSite | tuple[LatticeSite, ArrayLike]] | None = None,
                      planar_axis: ArrayLike | None = None,
                      field: ArrayLike | None = None,
                      step_size: float = 0.1,
                      initial_randomisation: InitialRandomisation | str = InitialRandomisation.JITTER,
                      seed: int | None = None,
                      rtol: float=1e-10,
                      atol: float=1e-12,
                      max_iters: int=1000,
                      verbose: bool=True):
        """ Get the classical ground state via gradient descent

        For more direct control, use the `ClassicalEnergyMinimisation` class

        :param fixed: List of sites that should be ignored by the minimisation, i.e. fixed in place
        :param planar: List of sites, or (site, axis) tuples that should be constrained to a plane
        :param planar_axis: Axis to use for planar constraints if not specified explicitly, default=[0,0,1]
        :param field: Magnetic field applied
        :param step_size: Size of step to make relative to the force, smaller values than the default might be needed
                          for systems with high energy.
        :param initial_randomisation: Option to RANDOMISE or JITTER the starting state, default JITTER, can also be NONE
        :param seed: Seed to use in any randomisation steps
        :param rtol: Convergence criterion, stop when [energy change] < rtol x [initial energy change]
        :param atol: Convergence criterion, stop when [energy change] < atol
        :param max_iters: Limit on the number of iterations
        :param verbose: Print information about the process to stdout

        :returns: A new `Hamiltonian` with optimised spin state
        """
        #
        # Build the constraints
        #

        # Deal with defaults
        default_axis = np.array([0,0,1]) if planar_axis is None else np.array(planar_axis)
        fixed = [] if fixed is None else fixed
        planar = [] if planar is None else planar

        # Mapping
        site_uid_lookup = {site.unique_id: index for index, site in enumerate(self.structure.sites)}


        # Actual building
        constraints = [Free for _ in self.structure.sites]

        for site in fixed:
            constraints[site_uid_lookup[site.unique_id]] = Fixed

        for site_or_tuple in planar:

            if isinstance(site_or_tuple, LatticeSite):
                index = site_uid_lookup[site_or_tuple.unique_id]
                constraints[index] = Planar(default_axis)

            elif isinstance(site_or_tuple, tuple):

                if len(site_or_tuple) != 2:
                    raise ValueError("Expected tuples in `planar` to be length 2")

                if not isinstance(site_or_tuple[0], LatticeSite):
                    raise TypeError("Expected first component of tuples to be a LatticeSite")

                try:
                    axis = np.array(site_or_tuple[1])
                except Exception:
                    raise TypeError("Expected second component of tuples to be an array")

                index = site_uid_lookup[site_or_tuple[0].unique_id]
                constraints[index] = Planar(axis)

            else:
                raise TypeError("Expected entries in planar to be sites, or tuples of sites and axes")

        # TODO: Apply constraints to symmetric sites?

        #
        # Run the minimisation
        #

        minimiser = ClassicalEnergyMinimisation(self, constraints, field, seed)
        minimiser.minimise(
            rtol=rtol,
            atol=atol,
            max_iters=max_iters,
            initial_randomisation=initial_randomisation,
            step_size=step_size,
            verbose=verbose)


        # Create new spins
        old_uid_to_new_site = {}
        new_sites = []
        for site_index, site in enumerate(minimiser.sites):
            spin_data = minimiser.moment_data[site_index, :, :]

            old_uid = site.unique_id
            new_site = LatticeSite(site.i, site.j, site.k,
                                   supercell_spins=spin_data,
                                   g=site.g, name=site.name)

            new_sites.append(new_site)
            old_uid_to_new_site[old_uid] = new_site

        structure = Structure(new_sites,
                              minimiser.hamiltonian.structure.unit_cell,
                              minimiser.hamiltonian.structure.spacegroup,
                              minimiser.hamiltonian.structure.supercell)

        exchanges = [exchange.updated(
            site_1=old_uid_to_new_site[exchange.site_1.unique_id],
            site_2=old_uid_to_new_site[exchange.site_2.unique_id])
            for exchange in minimiser.hamiltonian.exchanges]

        anisotropies = [anisotropy.updated(
            site=old_uid_to_new_site[anisotropy.site.unique_id])
            for anisotropy in minimiser.hamiltonian.anisotropies]

        return Hamiltonian(structure, exchanges, anisotropies)


    def _serialise(self, context: SPWSerialisationContext) -> dict:
        return {"magnetic_structure": self.structure._serialise(context),
                "exchanges": [exchange._serialise(context) for exchange in self.exchanges],
                "anisotropies": [anisotropy._serialise(context) for anisotropy in self.anisotopies]}

    @staticmethod
    @expects_keys("magnetic_structure, exchanges, anisotropies")
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        structure = Structure._deserialise(json["magnetic_structure"], context)
        exchanges = [Exchange._deserialise(exchange, context) for exchange in json["exchanges"]]
        anisotropies = [Anisotropy._deserialise(anisotropy, context) for anisotropy in json["anisotropies"]]

        return Hamiltonian(structure, exchanges, anisotropies)


class HamiltonianParameterization:
    """ Parameterisation of a Hamiltonian

    This is a callable class that returns a Hamiltonian with exchanges set by specified parameters
    """

    def __init__(self,
                 hamiltonian: Hamiltonian,
                 *parameters: ParametrizationType,
                 find_ground_state_with: dict | None = None):

        self._hamiltonian = hamiltonian
        self._find_ground_state = find_ground_state_with is not None
        self._ground_state_parameters = {} if find_ground_state_with is None else find_ground_state_with

        # Create a list of parameters that will be updated
        base_exchange_parameters: list[list[tuple[Exchange, str]]] = []
        base_anisotropy_parameters: list[list[tuple[Anisotropy, str]]] = []

        for param in parameters:
            exchanges, anisotropies = _regularise_parameters(hamiltonian, param)

            base_exchange_parameters.append(exchanges)
            base_anisotropy_parameters.append(anisotropies)

        # Check that there is no conflicts, this means that each parameter is only set by only one entry
        # Basically, we can just check for duplicates

        multicounts = [item for item, count in Counter(sum(base_exchange_parameters, [])).items() if count > 1]
        if len(multicounts) > 0:
            raise ValueError(f"Multiple parameters assigned to the same exchange parameter {multicounts}")

        multicounts = [item for item, count in Counter(sum(base_anisotropy_parameters, [])).items() if count > 1]
        if len(multicounts) > 0:
            raise ValueError(f"Multiple parameters assigned to the same anisotropy parameter {multicounts}")

        # Check that the exchanges have the required parameters
        for parameter_definition in base_exchange_parameters:
            for target, attribute in parameter_definition:
                if attribute not in target.parameters:
                    valid_parameters = ", ".join(target.parameters)
                    raise TypeError(f"{target} does not have parameter '{attribute}', it has: {valid_parameters}")


        # Check that the anisotropies have the required parameters
        for parameter_definition in base_anisotropy_parameters:
            for target, attribute in parameter_definition:
                if attribute not in target.scalar_parameters:
                    valid_parameters = ", ".join(target.scalar_parameters)
                    raise TypeError(f"{target} does not have parameter '{attribute}', it has: {valid_parameters}")

        # Convert the exchange to index in the list of exchanges
        exchange_unique_id_to_index = {exchange.unique_id: index
                                       for index, exchange in enumerate(hamiltonian.exchanges)}
        self._exchange_parameter_definitions: list[list[tuple[int, str]]] = []

        for parameter_definition in base_exchange_parameters:
            indexed_parameter_definition: list[tuple[int, str]] = []
            for target, attribute in parameter_definition:
                index = exchange_unique_id_to_index[target.unique_id]
                indexed_parameter_definition.append((index, attribute))
            self._exchange_parameter_definitions.append(indexed_parameter_definition)

        # Convert the anisotropy to index for the list of anisotropies
        anisotropy_unique_id_to_index = {anisotropy.unique_id: index
                                       for index, anisotropy in enumerate(hamiltonian.anisotropies)}
        self._anisotropy_parameter_definitions: list[list[tuple[int, str]]] = []

        for parameter_definition in base_anisotropy_parameters:
            indexed_parameter_definition: list[tuple[int, str]] = []
            for target, attribute in parameter_definition:
                index = anisotropy_unique_id_to_index[target.unique_id]
                indexed_parameter_definition.append((index, attribute))
            self._anisotropy_parameter_definitions.append(indexed_parameter_definition)

        # Number is useful
        self._n_parameters = len(self._exchange_parameter_definitions)

        assert len(self._anisotropy_parameter_definitions) == len(self._exchange_parameter_definitions), \
            "Should both be length n_parameters"

    @staticmethod
    def _legend_entry(exchange_name, parameter_name):
        if exchange_name is None or exchange_name == "":
            return parameter_name
        else:
            return exchange_name + "." + parameter_name

    def energy_plot(self,
                    parameter_values: ArrayLike,
                    path: Path,
                    field: ArrayLike | None = None,
                    show: bool = True,
                    colormap_name: str = 'jet',
                    new_figure: bool = True,
                    show_legend=True,
                    use_rust: bool = True):
        """ Show a plot of the energies """
        # Check/regularise input values
        parameter_values = np.array(parameter_values)
        if len(parameter_values.shape) == 1:
            parameter_values = parameter_values.reshape(-1, 1)

        if parameter_values.shape[1] != self._n_parameters:
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

        colors = cmap(np.linspace(0, 1, n_curves+2)[1:-1]) # Generally a bit nicer if we avoid the endpoints

        # Figure
        if new_figure:
            plt.figure("Energy")

        x_values = path.x_values()

        for i in range(n_curves):
            ham = self(*parameter_values[i, :])
            first = True
            for series in ham.sorted_positive_energies(path, field=field, use_rust=use_rust):
                if first:
                    label = ", ".join([f"{value:.3g}" for value in parameter_values[i, :]])
                    plt.plot(x_values, series, color=colors[i], label=label)
                    first=False
                else:
                    plt.plot(x_values, series, color=colors[i])

        if show_legend:
            plt.legend(loc="upper right")

        # Format the plot
        path.format_plot(plt)

        if show:
            plt.show()



    def __call__(self, *parameters: float) -> Hamiltonian:
        """ Get the Hamiltonian with parameters set """
        if len(parameters) != self._n_parameters:
            raise ValueError(f"Expected {self._n_parameters} parameters, got {len(parameters)}")

        # Updated exchanges
        new_exchanges = [exchange for exchange in self._hamiltonian.exchanges]
        for parameter_definition, value in zip(self._exchange_parameter_definitions, parameters):
            for (exchange_index, attribute) in parameter_definition:
                new_exchanges[exchange_index] = new_exchanges[exchange_index].updated(**{attribute: value})

        # Updated anisotropies
        new_anisotropies = [anisotropy for anisotropy in self._hamiltonian.anisotropies]
        for parameter_definition, value in zip(self._anisotropy_parameter_definitions, parameters):
            for anisotropy_index, attribute in parameter_definition:
                new_anisotropies[anisotropy_index] = new_anisotropies[anisotropy_index].updated(**{attribute: value})

        # Return new hamiltonian, optimise ground state if needed
        new_hamiltonian = Hamiltonian(self._hamiltonian.structure, new_exchanges, self._hamiltonian.anisotropies)
        if self._find_ground_state:
            return new_hamiltonian.ground_state(**self._ground_state_parameters)
        else:
            return new_hamiltonian

    def __repr__(self):
        parts = []
        for index, (exchanges, anisotropies) in enumerate(zip(self._exchange_parameter_definitions,
                                                              self._anisotropy_parameter_definitions)):

            exchange_parts = [f"{self._hamiltonian.exchanges[index].name}.{parameter}"
                                for index, parameter in exchanges]

            anisotropy_parts = [f"{self._hamiltonian.anisotropies[index]}.{parameter}"
                                for index, parameter in anisotropies]

            data = ", ".join(exchange_parts + anisotropy_parts)

            parts.append(f"argument_{index} -> {data}")
        s = "; ".join(parts)

        return f"HamiltonianParameterization({s})"
