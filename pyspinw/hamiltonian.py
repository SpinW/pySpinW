"""Selection of different Hamiltonians"""
import logging
from abc import ABC, abstractmethod
from collections.abc import Callable

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
from pyspinw.symmetry.supercell import TrivialSupercell
from pyspinw.units import IntensityUnits, intensity_units


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
    # Ensure there are the same number of modes throughout
    nans = np.isnan(en_out)
    nanmodes = np.sum(nans, axis=0)
    n_modes = len(np.where((nanmodes < en_out.shape[0]))[0])
    for col in set([int(idx) for accidentals in np.where((nanmodes > 0) * (nanmodes < en_out.shape[0]))[0]
                    for idx in np.where(nans[:,accidentals])[0]]):
        idnxt = 1 if col < en_out.shape[0]-1 else -1
        idnan = np.where(~np.isnan(en_out[col + idnxt,:]))[0]
        idx = np.array([np.nanargmin(np.abs(en_out[col + idnxt,i] - en_out[col]))
               for i in idnan])
        en_out[col,idnan] = np.take(en_out[col, :], idx)
        int_out[col,idnan] = np.take(int_out[col, np.where(~np.isnan(int_out[col,:]))] / np.bincount(idx), idx)
    # Removes columns which are all NaNs
    en_out = np.delete(en_out, np.where(nanmodes == en_out.shape[0])[0], axis=1)
    int_out = np.delete(int_out, np.where(nanmodes == en_out.shape[0])[0], axis=1)
    return en_out, int_out


def egrid(energy: ArrayLike, intensity: ArrayLike, evect: ArrayLike, dE: float | Callable):
    """Bins a set of energy/intensity into a spectrum with Gaussian broadening"""
    esigma = dE(evect) if isinstance(dE, Callable) else np.ones(evect.shape) * dE
    esigma /= 2 * np.sqrt(2* np.log(2))  # Convert from FWHM
    spec = np.zeros((energy.shape[0], len(evect)))
    ew_norm = np.array(list(np.diff(evect)) + [evect[-1] - evect[-2]]) / (np.sqrt(2 * np.pi) * esigma)
    for iE in range(energy.shape[1]):
        spec += intensity[:,iE,np.newaxis] * ew_norm * np.exp(-0.5 * ((evect[np.newaxis,:] - energy[:,iE,np.newaxis]) / esigma)**2)
    spec[np.where(spec < np.finfo(np.float32).eps)] = np.nan
    return spec


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
        if intensity_unit == IntensityUnits.PERCELL:
            scale_factor = rotations.shape[0] / np.prod(scaling)
            intensity = [res * scale_factor for res in result[1]]

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
