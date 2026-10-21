"""Different kinds of samples."""

from abc import ABC, abstractmethod
from enum import Enum
from typing import Sequence, Callable

import numpy as np
from numpy._typing import ArrayLike

import matplotlib.pyplot as plt

from pyspinw.serialisation import SPWSerialisable, SPWSerialisationContext, SPWDeserialisationContext, expects_keys, \
    numpy_serialise, numpy_deserialise
from pyspinw.units import IntensityUnits, CoordsUnits
from pyspinw.calculations.spherical_integration import SphericalPointGeneratorType, point_generator
from pyspinw.checks import check_sizes
from pyspinw.hamiltonian import Hamiltonian, ParametrizationType, omegasum, egrid
from pyspinw.path import Path, Path1D, Path1DBase, EmpiricalPath1D, Slice
from pyspinw.tolerances import tolerances


# pylint: disable=R0903

class Sample(ABC, SPWSerialisable):
    """Representation of the macrostructure of a sample used in an experiment (Twin, Powder etc).

    Parameters
    ----------
    hamiltonian : Hamiltonian
        Spin Hamiltonian of the system.

    Attributes
    ----------
    hamiltonian : Hamiltonian
        Spin Hamiltonian of the system.
    """

    serialisation_name = "sample"
    sample_name = "<base-class>"

    def __init__(self, hamiltonian: Hamiltonian):
        self.hamiltonian = hamiltonian

    def _energies(self,
                 path: Path,
                 field: ArrayLike | None = None,
                 use_rust: bool=True):
        """Calculate sorted positive magnon energies along a path.

        Parameters
        ----------
        path : Path
            Reciprocal-space path along which to evaluate energies.
        field : array_like, optional
            External magnetic field vector of shape ``(3,)``.
        use_rust : bool, default=True
            Whether to use the Rust backend for calculation.

        Returns
        -------
        list of numpy.ndarray
            Sorted positive mode energies at each q-point.
        """
        return self.hamiltonian.sorted_positive_energies(path, field=field, use_rust=use_rust)

    def _serialise(self, context: SPWSerialisationContext) -> dict:
        return {"type": self.sample_name,
                "data": self._sample_serialise(context)}

    @staticmethod
    @expects_keys("type,data")
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        return _all_samples_lookup[json["type"]]._sample_deserialise(json["data"], context)

    @abstractmethod
    def _sample_serialise(self, context: SPWSerialisationContext):
        """Serialise the data for specific kind of sample.

        Parameters
        ----------
        context : SPWSerialisationContext
            Serialisation context.

        Returns
        -------
        dict
            Serialised sample data.
        """

    @staticmethod
    @abstractmethod
    def _sample_deserialise(json, context: SPWDeserialisationContext):
        """Deserialise the data for this kind of sample.

        Parameters
        ----------
        json : dict
            Serialised sample dictionary.
        context : SPWDeserialisationContext
            Deserialisation context.

        Returns
        -------
        Sample
            Deserialised sample instance.
        """

    def _serialise(self, context: SPWSerialisationContext) -> dict:
        return {"type": self.sample_name,
                "data": self._sample_serialise(context)}

    @staticmethod
    @expects_keys("type,data")
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        return _all_samples_lookup[json["type"]]._sample_deserialise(json["data"], context)

    @abstractmethod
    def _sample_serialise(self, context: SPWSerialisationContext):
        """Serialise the data for specific kind of sample.

        Parameters
        ----------
        context : SPWSerialisationContext
            Serialisation context.

        Returns
        -------
        dict
            Serialised sample data.
        """

    @staticmethod
    @abstractmethod
    def _sample_deserialise(json, context: SPWDeserialisationContext):
        """Deserialise the data for this kind of sample.

        Parameters
        ----------
        json : dict
            Serialised sample dictionary.
        context : SPWDeserialisationContext
            Deserialisation context.

        Returns
        -------
        Sample
            Deserialised sample instance.
        """

class Sample3D(Sample):
    """Sample where the direction of q matters."""

    @abstractmethod
    def _energies_and_intensities(self,
                                 points: np.ndarray,
                                 field: ArrayLike | None = None,
                                 use_rust: bool = True,
                                 use_rotating: bool=True,
                                 intensity_unit: IntensityUnits | str='cell',):
        """Abstract method to get the energies and intensities.

        Parameters
        ----------
        points : numpy.ndarray
            Array of 3D q-points of shape ``(n_q, 3)``.
        field : array_like, optional
            External magnetic field vector of shape ``(3,)``.
        use_rust : bool, default=True
            Whether to use the Rust backend for calculation.
        use_rotating : bool, default=True
            Whether to use rotating-frame coordinates.
        intensity_unit : IntensityUnits or str, default='cell'
            Unit system for calculated intensities.

        Returns
        -------
        energies : numpy.ndarray
            Calculated mode energies.
        intensities : list or numpy.ndarray
            Calculated mode intensities.
        """

    @check_sizes(field=(3,), allow_nones=True)
    def energies(self,
                 path: Path,
                 field: ArrayLike | None = None,
                 use_rust: bool = True) -> list[np.ndarray]:
        """Get the energies along a specified path.

        Parameters
        ----------
        path : Path
            Reciprocal-space path along which to evaluate energies.
        field : array_like, optional
            External magnetic field vector of shape ``(3,)``.
        use_rust : bool, default=True
            Whether to use the Rust backend for calculation.

        Returns
        -------
        list of numpy.ndarray
            Sorted positive mode energies along `path`.
        """
        return self._energies(path, field, use_rust)


    def energies_and_intensities(self,
                                 path: Path,
                                 field: ArrayLike | None = None,
                                 use_rust: bool = True,
                                 use_rotating: bool = True,
                                 intensity_unit: IntensityUnits | str = 'cell'):
        """Get energy and intensity data.

        Parameters
        ----------
        path : Path
            Reciprocal-space path along which to evaluate energies and intensities.
        field : array_like, optional
            External magnetic field vector of shape ``(3,)``.
        use_rust : bool, default=True
            Whether to use the Rust backend for calculation.
        use_rotating : bool, default=True
            Whether to use rotating-frame coordinates.
        intensity_unit : IntensityUnits or str, default='cell'
            Unit system for calculated intensities.

        Returns
        -------
        energies : numpy.ndarray
            Calculated mode energies at each q-point along `path`.
        intensities : list or numpy.ndarray
            Calculated mode intensities at each q-point along `path`.
        """
        return self._energies_and_intensities(
            path.q_points(),
            field=field,
            use_rust=use_rust,
            use_rotating=use_rotating,
            intensity_unit=intensity_unit)

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
                       intensity_unit: IntensityUnits | str = 'cell'):
        """Create a spaghetti diagram with intensity as colorfill overplotted by mode energies.

        Parameters
        ----------
        path : Path
            Reciprocal-space path to plot.
        evect : array_like, optional
            1D array of energy bin edges or centres in meV.
        dE : float or callable, optional
            Energy resolution broadening (FWHM in meV) or resolution function.
        vmin : float, default=0
            Minimum value for the color scale.
        vmax : float, optional
            Maximum value for the color scale.
        field : array_like, optional
            External magnetic field vector of shape ``(3,)``.
        show : bool, default=True
            If ``True``, call ``matplotlib.pyplot.show()``.
        new_figure : bool, default=True
            If ``True``, create a new matplotlib figure and axes.
        use_rust : bool, default=True
            Whether to use the Rust backend for calculation.
        use_rotating : bool, default=True
            Whether to use rotating-frame coordinates.
        intensity_unit : IntensityUnits or str, default='cell'
            Unit system for calculated intensities.
        scale : str, default='linear'
            Intensity scaling mode.

        Returns
        -------
        matplotlib.figure.Figure or None
            The matplotlib figure if ``show=False``, otherwise ``None``.
        """
        if new_figure:
            fig, ax = plt.subplots()
        else:
            fig = plt.gcf()
            ax = fig.get_axes()[0]

        x_values = path.x_values()
        energy, intensity = omegasum(*self._energies_and_intensities(path.q_points(), field,
                                                                    use_rust, use_rotating, intensity_unit))
        if evect is None:
            emax = np.max([0.001, np.nanmax(np.real(energy))])
            denom = 10**np.floor(np.log10(emax))
            evect = np.linspace(0, np.ceil(emax / denom) * denom, 100)
        if dE is None:
            dE = np.mean(np.diff(evect)) * 2.35
        spec = egrid(energy, intensity, evect, dE)

        if vmax is None:
            try:
                vmax = np.nanmax(spec[:, np.where(evect > max(np.max(evect)/10, 0.1))]) / 10.
            except ValueError:
                vmax = 1.0

        if dE != 0:
            mesh = ax.pcolormesh(x_values, evect, spec.T, vmin=vmin, vmax=vmax, cmap="magma_r")
            fig.colorbar(mesh, ax=ax)

        ax.plot(x_values, np.real(energy), 'k')
        if (np.imag(energy) > 0.01).any():
            ax.plot(x_values, np.imag(energy), 'or')
        ax.set_ylim(0, np.max(evect))
        path.format_plot(ax)
        ax.set_ylabel('Magnon Energy (meV)')

        if show:
            plt.show()
        else:
            return fig

    def intensity_map(self,
                            q_slice: Slice,
                            vmin: float = 0,
                            vmax: float | None = None,
                            field: ArrayLike | None = None,
                            show: bool = True,
                            new_figure: bool = True,
                            use_rust: bool = True,
                            use_rotating: bool = True,
                            intensity_unit: IntensityUnits | str = 'cell',
                            cmap: str = "magma_r",
                            origin: str = "lower"):
        """Plot a reciprocal space map with intensity integrated over energy.

        Parameters
        ----------
        q_slice : Slice
            2D slice object defining the grid of q-points. Energy bounds
            ``e_min`` and ``e_max`` (in meV) are taken from this object.
        vmin : float, default=0
            Minimum limit of the color scale.
        vmax : float, optional
            Maximum limit of the color scale.
        field : array_like, optional
            External magnetic field vector of shape ``(3,)``.
        show : bool, default=True
            Whether to display the plot.
        new_figure : bool, default=True
            Whether to create a new figure.
        use_rust : bool, default=True
            Whether to use the Rust backend for calculation.
        use_rotating : bool, default=True
            Whether to use rotating-frame coordinates.
        intensity_unit : IntensityUnits or str, default='cell'
            Unit system for calculated intensities.
        cmap : str, default='magma_r'
            Matplotlib colormap name.
        origin : str, default='lower'
            ``'lower'`` or ``'upper'`` for matplotlib ``imshow`` origin.

        Returns
        -------
        matplotlib.figure.Figure or None
            The matplotlib figure if ``show=False``, otherwise ``None``.
        """
        if new_figure:
            fig, ax = plt.subplots()
        else:
            fig = plt.gcf()
            ax = fig.get_axes()[0]

        # Get energy bounds from the Slice object
        e_min = q_slice.e_min
        e_max = q_slice.e_max

        # Step 1: Calculate all energies and intensities at each q-point
        energies, intensities = self._energies_and_intensities(
            q_slice.q_points(),
            field=field,
            use_rust=use_rust,
            use_rotating=use_rotating,
            intensity_unit=intensity_unit
        )

        # Step 2: Convert to numpy arrays and nan is counted as 0
        energies = np.array(energies)
        intensities = np.array(intensities)

        energies = np.nan_to_num(energies, nan=0.0)
        intensities = np.nan_to_num(intensities, nan=0.0)


        # Step 3: Energy integration - mask and sum
        mask = energies > 0  # Only positive energy modes
        if e_min is not None:
            mask &= energies >= e_min
        if e_max is not None:
            mask &= energies <= e_max

        # Sum over all modes within the energy window
        integrated = np.sum(intensities * mask, axis=1)  # shape = (n_q,)

        # Step 4: Reshape back to 2D grid
        # meshgrid(indexing="xy") lays out data as (n_b, n_a) in memory
        # imshow expects shape (rows, cols) = (n_b, n_a) where rows = vertical/axis_2
        grid_data = integrated.reshape(q_slice.n_b, q_slice.n_a)

        # Step 5: Auto-scale color limits if not provided
        if vmax is None:
            vmax = np.nanmax(grid_data)

        # Step 6: Plot heatmap
        im = ax.imshow(
            grid_data,
            extent=q_slice.extent(),
            origin=origin,
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
            aspect="auto"
        )

        # Step 7: Add colorbar and labels
        fig.colorbar(im, ax=ax, label="Intensity")
        q_slice.format_plot(ax)

        # Add title showing energy window
        if e_min is not None and e_max is not None:
            ax.set_title(f"Reciprocal Space Map (E = {e_min:.0f}-{e_max:.0f} meV)")
        elif e_min is not None:
            ax.set_title(f"Reciprocal Space Map (E > {e_min:.0f} meV)")
        elif e_max is not None:
            ax.set_title(f"Reciprocal Space Map (E < {e_max:.0f} meV)")
        else:
            ax.set_title("Reciprocal Space Map (all energies)")

        fig.tight_layout()

        if show:
            plt.show()
        else:
            return fig


class Sample1D(Sample):
    """Sample where the direction of q does not matter.

    An important difference from Sample3D is there's no option to just get energies,
    this is because they will form a continuum in 1D. The 3D q-energy curves can
    be obtained from the Hamiltonian, they're just not exposed at this level.
    """



class CrystalDomain(SPWSerialisable):
    """Orientation and amount of crystalline domain (e.g. one of a twin, or a domain in the traditional sense).

    Parameters
    ----------
    transformation : (3, 3) array_like
        Orthonormal ``3x3`` transformation matrix describing the domain orientation.
    weighting : float
        Non-negative weight (relative volume fraction) of the domain.

    Attributes
    ----------
    transformation : (3, 3) numpy.ndarray
        Orthonormal ``3x3`` transformation matrix describing the domain orientation.
    weighting : float
        Non-negative weight of the domain.

    Raises
    ------
    ValueError
        If `transformation` is not orthonormal or if `weighting` is negative.
    """

    serialisation_name = "crystal_domain"

    @check_sizes(transformation=(3,3), force_numpy=True)
    def __init__(self, transformation: ArrayLike, weighting: float):

        transformation = np.array(transformation, dtype=float)

        # check the transformation is valid
        if not np.allclose(transformation @ transformation.T, np.eye(3), atol=tolerances.IS_ZERO_TOL):
            raise ValueError(f"Expected crystal domain transformation to be orthonormal, got {transformation}")

        if weighting < 0:
            raise ValueError(f"Weightings must be positive, got {weighting}")

        self.transformation = transformation
        self.weighting = weighting

    def _serialise(self, context: SPWSerialisationContext):
        return {
            "transformation": numpy_serialise(self.transformation),
            "weighting": float(self.weighting)
        }

    @staticmethod
    @expects_keys("transformation,weighting")
    def _deserialise(json: dict, context: SPWDeserialisationContext):
        return CrystalDomain(
            transformation=numpy_deserialise(json["transformation"]),
            weighting=json["weighting"])

CrystalDomainLike = CrystalDomain | tuple[ArrayLike, float]


def _domain_like_to_domain(domain_like: CrystalDomainLike) -> CrystalDomain:
    """Convert CrystalDomainLike to CrystalDomain.

    Parameters
    ----------
    domain_like : CrystalDomain or tuple of (array_like, float)
        Domain object or ``(transformation, weighting)`` tuple.

    Returns
    -------
    CrystalDomain
        Corresponding :class:`CrystalDomain` instance.
    """
    if isinstance(domain_like, CrystalDomain):
        return domain_like
    else:
        return CrystalDomain(*domain_like)


class SingleCrystal(Sample3D):
    """Specifies a single crystal sample.

    Parameters
    ----------
    hamiltonian : Hamiltonian
        Spin Hamiltonian of the single crystal.
    """

    sample_name = "single_crystal"

    def __init__(self, hamiltonian: Hamiltonian):

        super().__init__(hamiltonian=hamiltonian)

    def _energies_and_intensities(self,
                                 q_points: np.ndarray,
                                 field: ArrayLike | None = None,
                                 use_rust: bool = True,
                                 use_rotating: bool = True,
                                 intensity_unit: IntensityUnits | str = 'cell'):
        """Get energy and intensity data.

        Parameters
        ----------
        q_points : numpy.ndarray
            Array of 3D q-points of shape ``(n_q, 3)``.
        field : array_like, optional
            External magnetic field vector of shape ``(3,)``.
        use_rust : bool, default=True
            Whether to use the Rust backend for calculation.
        use_rotating : bool, default=True
            Whether to use rotating-frame coordinates.
        intensity_unit : IntensityUnits or str, default='cell'
            Unit system for calculated intensities.

        Returns
        -------
        energies : numpy.ndarray
            Calculated mode energies.
        intensities : list or numpy.ndarray
            Calculated mode intensities.
        """
        # Note: we don't use the _ method because we want to ignore non-magnetic sites
        return self.hamiltonian._energies_and_intensities(
            q_points, field=field, use_rust=use_rust, use_rotating=use_rotating, intensity_unit=intensity_unit)

    def _sample_serialise(self, context: SPWSerialisationContext):
        return {"hamiltonian": self.hamiltonian._serialise(context)}

    @staticmethod
    @expects_keys("hamiltonian")
    def _sample_deserialise(json, context: SPWDeserialisationContext):
        return SingleCrystal(Hamiltonian._deserialise(json["hamiltonian"], context))


class Multidomain(Sample3D):
    """Sample consisting of multiple domains.

    Parameters
    ----------
    hamiltonian : Hamiltonian
        Spin Hamiltonian of the system.
    domains : list of CrystalDomainLike
        List of crystalline domains or ``(transformation, weighting)`` tuples.

    Attributes
    ----------
    weights : list of float
        Normalised domain weights summing to ``1.0``.
    """

    sample_name = "multidomain"

    def __init__(self,
                 hamiltonian: Hamiltonian,
                 domains: list[CrystalDomainLike]):

        self._domains = [_domain_like_to_domain(domain) for domain in domains]

        self._transformations = [domain.transformation for domain in self._domains]

        weights = [domain.weighting for domain in self._domains]
        total_weight = sum(weights)
        self.weights = [weight/total_weight for weight in weights]

        super().__init__(hamiltonian=hamiltonian)


    def _energies(self,
                 path: Path,
                 field: ArrayLike | None = None,
                 use_rust: bool=True) -> list[np.ndarray]:

        q_vectors = path.q_points()

        output = []
        for transformation in self._transformations:

            transformed_q = q_vectors @ transformation.T
            transformed_field = None if field is None else transformation @ field

            output += self.hamiltonian._sorted_positive_energies(
                        transformed_q,
                        field=transformed_field,
                        use_rust=use_rust)

        return output

    def _energies_and_intensities(self,
                                  q_vectors: np.ndarray,
                                  field: ArrayLike | None = None,
                                  use_rust: bool = True,
                                  use_rotating: bool = True,
                                  intensity_unit: IntensityUnits | str = 'cell'):

        n_q = q_vectors.shape[0]

        if field is not None:
            field = np.array(field).reshape(3)

        output_energies = []
        output_intensities = [[] for _ in range(n_q)]
        for transformation, weight in zip(self._transformations, self.weights):

            transformed_q = q_vectors @ transformation.T
            transformed_field = None if field is None else field @ transformation.T

            energy, intensities = self.hamiltonian._energies_and_intensities(
                        transformed_q,
                        field=transformed_field,
                        use_rust=use_rust)

            output_energies.append(energy)
            for i in range(n_q):
                output_intensities[i] += [weight * intensity for intensity in intensities[i]]

        output_energies = np.concatenate(output_energies, axis=1)

        return output_energies, output_intensities

    @staticmethod
    @expects_keys("hamiltonian,domains")
    def _sample_deserialise(json, context: SPWDeserialisationContext):
        return Multidomain(
            hamiltonian=Hamiltonian._deserialise(json["hamiltonian"], context),
            domains=[CrystalDomain._deserialise(data, context) for data in json["domains"]]
        )

    def _sample_serialise(self, context: SPWSerialisationContext):
        return {
            "hamiltonian": self.hamiltonian._serialise(context),
            "domains": [domain._serialise(context) for domain in self._domains]
        }


class Twin(Multidomain):
    """Specify a twinned crystal by a twinning plane.

    This is a special case of :class:`Multidomain`, which can use completely general transformations,
    and have more than one "twin".

    Parameters
    ----------
    hamiltonian : Hamiltonian
        Spin Hamiltonian of the system.
    twinning_plane_normal : (3,) array_like
        Normal vector to the twinning mirror plane.
    second_twin_fraction : float
        Volume fraction of the second twin domain, between ``0`` and ``1``.
    units : CoordsUnits or str, default='lu'
        Coordinate units of `twinning_plane_normal` (``'lu'`` or ``'xyz'``).

    Raises
    ------
    ValueError
        If `second_twin_fraction` is not in ``[0, 1]`` or if `units` is not
        recognised.
    """

    sample_name = "twin"

    @check_sizes(twinning_plane_normal=(3,), force_numpy=True)
    def __init__(self,
                 hamiltonian: Hamiltonian,
                 twinning_plane_normal: ArrayLike,
                 second_twin_fraction: float,
                 units: CoordsUnits | str = "lu"):

        if not (0 <= second_twin_fraction <= 1):
            raise ValueError("Expected weighting to be between 0 and 1")

        self._twinning_plane_normal = twinning_plane_normal
        self._second_twin_fraction = second_twin_fraction

        if isinstance(units, str):
            units = CoordsUnits(units)

        self._units = units

        # Get the transformation in lattice units
        match units:
            case CoordsUnits.LU:
                vector = twinning_plane_normal
            case CoordsUnits.XYZ:
                vector = hamiltonian.structure.unit_cell.cartesian_to_lattice_units(twinning_plane_normal)
            case _:
                raise ValueError(f"Unknown unit type {units}")

        # Vector needs to be normalised
        vector = np.array(vector, dtype=float)
        vector /= np.sqrt(np.sum(vector**2))

        # Reflection in plane specified by the normal
        transform = np.eye(3) - 2 * vector.reshape(1, 3) * vector.reshape(3, 1)

        first_twin_fraction = 1 - second_twin_fraction

        super().__init__(
            hamiltonian=hamiltonian,
            domains=[CrystalDomain(np.eye(3), first_twin_fraction),
                     CrystalDomain(transform, second_twin_fraction)])

    def _sample_serialise(self, context: SPWSerialisationContext):
        return {
            "hamiltonian": self.hamiltonian._serialise(context),
            "normal": numpy_serialise(self._twinning_plane_normal),
            "fraction": self._second_twin_fraction,
            "units": self._units.value
        }

    @staticmethod
    @expects_keys("hamiltonian,normal,fraction,units")
    def _sample_deserialise(json, context: SPWDeserialisationContext):
        return Twin(
            hamiltonian=Hamiltonian._deserialise(json["hamiltonian"], context),
            twinning_plane_normal=numpy_deserialise(json["normal"]),
            second_twin_fraction=json["fraction"],
            units=CoordsUnits(json["units"])
        )

class ScalingMethod(Enum):
    """Scaling methods for plots."""

    LINEAR = 'linear'
    LOG = 'log'

class ParameterizedPowderSpectrum:
    """A powder spectrum that has been parameterised.

    Calling an instance with parameter values evaluates the powder spectrum at
    those specified parameter values.

    Parameters
    ----------
    powder : Powder
        Powder sample to parameterise.
    parameters : Sequence of ParametrizationType
        Parameter specifications passed to :meth:`Hamiltonian.parameterize`.
    path : Path1D or EmpiricalPath1D or array_like
        1D path or array of ``|q|`` magnitudes in inverse Angstroms.
    n_samples : int, default=100
        Approximate number of spherical orientation samples per ``|q|`` shell.
    method : SphericalPointGeneratorType or str, default='fibonacci'
        Spherical point sampling method.
    min_energy : float, optional
        Minimum energy for histogramming in meV.
    max_energy : float, optional
        Maximum energy for histogramming in meV.
    n_energy_bins : int, optional
        Number of energy bins.
    energy_stddev : float, optional
        Standard deviation for Gaussian energy broadening in meV.
    random_seed : int, optional
        Random seed for stochastic spherical sampling methods.
    use_rust : bool, default=True
        Whether to use the Rust backend for calculation.
    find_ground_state_with : dict, optional
        Keyword arguments forwarded to ground-state energy minimisation when
        updating parameters.
    ignore_imaginary : bool, default=False
        If ``True``, zero out intensities of modes with significant imaginary
        energy components.

    Attributes
    ----------
    parameterised_hamiltonian : callable
        Parameterized Hamiltonian factory returned by
        :meth:`Hamiltonian.parameterize`.
    """

    def __init__(self,
                 powder: "Powder",
                 parameters: Sequence[ParametrizationType],
                 path: Path1D | EmpiricalPath1D | ArrayLike,
                 n_samples: int = 100,
                 method: SphericalPointGeneratorType | str = "fibonacci",
                 min_energy: float | None = None,
                 max_energy: float | None = None,
                 n_energy_bins: int | None = None,
                 energy_stddev: float | None = None,
                 random_seed: int | None = None,
                 use_rust: bool = True,
                 find_ground_state_with: dict | None=None,
                 ignore_imaginary: bool = False):

        self._path = path
        self._n_samples = n_samples
        self._method = method
        self._min_energy = min_energy
        self._max_energy = max_energy
        self._n_energy_bins = n_energy_bins
        self._energy_stddev = energy_stddev
        self._random_seed = random_seed
        self._use_rust = use_rust
        self._find_ground_state_with = find_ground_state_with
        self._ignore_imaginary = ignore_imaginary

        self.parameterised_hamiltonian = powder.hamiltonian.parameterize(
                                *parameters,
                                           find_ground_state_with=find_ground_state_with)

    def __call__(self, *parameters):  # noqa: D102
        powder = Powder(self.parameterised_hamiltonian(*parameters))
        return powder.spectrum(
            path=self._path,
            n_samples=self._n_samples,
            method=self._method,
            min_energy=self._min_energy,
            max_energy=self._max_energy,
            n_energy_bins=self._n_energy_bins,
            energy_stddev=self._energy_stddev,
            random_seed=self._random_seed,
            use_rust=self._use_rust,
            ignore_imaginary=self._ignore_imaginary)[2]


class Powder(Sample1D):
    """Sample is a powder.

    Parameters
    ----------
    hamiltonian : Hamiltonian
        Spin Hamiltonian of the system.
    """

    def __init__(self, hamiltonian: Hamiltonian):

        super().__init__(hamiltonian)

    def spectrum(self,
                 path: Path1D | EmpiricalPath1D | ArrayLike,
                 n_samples: int=100,
                 method: SphericalPointGeneratorType | str = "fibonacci",
                 min_energy: float | None = None,
                 max_energy: float | None = None,
                 n_energy_bins: int | None = None,
                 energy_stddev: float | None = None,
                 random_seed: int | None = None,
                 use_rust: bool = True,
                 ignore_imaginary: bool = False):
        """Get the powder spectrum.

        Parameters
        ----------
        path : Path1D or EmpiricalPath1D or array_like
            1D path or array of ``|q|`` magnitudes in inverse Angstroms.
        n_samples : int, default=100
            Approximate number of spherical orientation samples per ``|q|`` shell.
        method : SphericalPointGeneratorType or str, default='fibonacci'
            Spherical point sampling method.
        min_energy : float, optional
            Minimum energy for histogramming in meV.
        max_energy : float, optional
            Maximum energy for histogramming in meV.
        n_energy_bins : int, optional
            Number of energy bins. Defaults to ``n_samples // 4``.
        energy_stddev : float, optional
            Standard deviation for Gaussian energy broadening in meV. If ``None``,
            histogram binning is used.
        random_seed : int, optional
            Random seed for stochastic spherical sampling methods.
        use_rust : bool, default=True
            Whether to use the Rust backend for calculation.
        ignore_imaginary : bool, default=False
            If ``True``, zero out intensities of modes with significant imaginary
            energy components.

        Returns
        -------
        q_values : numpy.ndarray
            1D array of ``|q|`` magnitudes.
        energy_bin_centres : numpy.ndarray
            1D array of energy bin centres in meV.
        output : numpy.ndarray
            2D array of powder-averaged intensities of shape
            ``(len(q_values), n_energy_bins)``.
        """
        if not isinstance(path, Path1DBase):
            path = EmpiricalPath1D(path)

        generator = point_generator(method)(n_samples, seed = random_seed)

        chunk_size = generator.actual_n_points
        qs = path.q_values()

        # The sphere has radius |q| in inverse angstrom, which is what a powder measures, but the
        # calculation takes lattice units: q_lu = q_xyz A^T / 2 pi, where the rows of A are the
        # real-space lattice vectors
        cartesian_to_lattice_units = self.hamiltonian.structure.unit_cell._xyz.T / (2 * np.pi)
        points = np.concatenate([(generator.points * q) @ cartesian_to_lattice_units for q in qs], axis=0)

        energies, intensities = self.hamiltonian._energies_and_intensities(points, use_rust=use_rust)

        intensities = np.real(np.array(intensities))
        if ignore_imaginary:
            i_imag = np.where(np.abs(np.imag(energies)) > tolerances.IMAG_MODE_TOL)
            intensities[i_imag] = 0
        energies = np.real(np.array(energies))

        positive_energies = energies > 0


        # Historgramming parameters
        min_energy = np.min(energies[positive_energies]) if min_energy is None else min_energy
        max_energy = np.max(energies[positive_energies]) if max_energy is None else max_energy
        n_energy_bins = n_samples // 4 if n_energy_bins is None else n_energy_bins

        energy_bin_edges = np.linspace(min_energy, max_energy, n_energy_bins + 1)
        energy_bin_centres = 0.5 * (energy_bin_edges[1:] + energy_bin_edges[:-1])

        # Output map
        output = np.zeros((len(qs), n_energy_bins))

        if energy_stddev is not None:
            # Gausian scaling
            exponent_factor = -0.5 / (energy_stddev**2)
            # Gaussian normalisation
            normalisation_factor = 1.0 / np.sqrt(2*np.pi*(energy_stddev**2))

        # at each q, bin the energies with weights of the intensities, but ignore the negative ones
        for i, q in enumerate(path.q_values()):

            # Get the part of the data that's relevant
            energy = energies[i*chunk_size:(i+1)*chunk_size]
            intensity = intensities[i*chunk_size:(i+1)*chunk_size]
            positives = positive_energies[i*chunk_size:(i+1)*chunk_size]

            energy = energy[positives]
            intensity = intensity[positives]


            if energy_stddev is None:
                # Use a binning method
                binned, bin_edges = np.histogram(energy, bins=energy_bin_edges, weights=intensity)
            else:
                binned = np.zeros((n_energy_bins,))
                # Add up gaussians
                for energy_value, intensity_value in zip(energy, intensity):
                    gaussian = np.exp(((energy_bin_centres - energy_value)**2)*exponent_factor)
                    binned += (intensity_value * normalisation_factor) * gaussian


            output[i, :] = binned / chunk_size

        return path.q_values(), energy_bin_centres, output

    def parameterized_spectrum(self,
            parameters: Sequence[ParametrizationType],
            path: Path1D | EmpiricalPath1D | ArrayLike,
            n_samples: int = 100,
            method: SphericalPointGeneratorType | str = "fibonacci",
            min_energy: float | None = None,
            max_energy: float | None = None,
            n_energy_bins: int | None = None,
            energy_stddev: float | None = None,
            random_seed: int | None = None,
            use_rust: bool = True,
            find_ground_state_with: dict | None = None,
            ignore_imaginary: bool = False):
        """Get the powder spectrum as a function of parameters.

        Parameters
        ----------
        parameters : Sequence of ParametrizationType
            Parameter specifications passed to :meth:`Hamiltonian.parameterize`.
        path : Path1D or EmpiricalPath1D or array_like
            1D path or array of ``|q|`` magnitudes in inverse Angstroms.
        n_samples : int, default=100
            Approximate number of spherical orientation samples per ``|q|`` shell.
        method : SphericalPointGeneratorType or str, default='fibonacci'
            Spherical point sampling method.
        min_energy : float, optional
            Minimum energy for histogramming in meV.
        max_energy : float, optional
            Maximum energy for histogramming in meV.
        n_energy_bins : int, optional
            Number of energy bins.
        energy_stddev : float, optional
            Standard deviation for Gaussian energy broadening in meV.
        random_seed : int, optional
            Random seed for stochastic spherical sampling methods.
        use_rust : bool, default=True
            Whether to use the Rust backend for calculation.
        find_ground_state_with : dict, optional
            Keyword arguments forwarded to ground-state energy minimisation when
            updating parameters.
        ignore_imaginary : bool, default=False
            If ``True``, zero out intensities of modes with significant imaginary
            energy components.

        Returns
        -------
        ParameterizedPowderSpectrum
            Callable parameterized powder spectrum object.
        """
        return ParameterizedPowderSpectrum(
                        self,
                        parameters,
                        path=path,
                        n_samples=n_samples,
                        method=method,
                        min_energy=min_energy,
                        max_energy=max_energy,
                        n_energy_bins=n_energy_bins,
                        energy_stddev=energy_stddev,
                        random_seed=random_seed,
                        use_rust=use_rust,
                        find_ground_state_with=find_ground_state_with,
                        ignore_imaginary=ignore_imaginary)

    def show_spectrum(self,
                      path: Path1D | EmpiricalPath1D | ArrayLike | tuple,
                      n_samples: int = 100,
                      method: SphericalPointGeneratorType | str = "fibonacci",
                      min_energy: float | None = None,
                      max_energy: float | None = None,
                      n_energy_bins: int | None = None,
                      energy_stddev: float | None = None,
                      random_seed: int | None = None,
                      scaling_method: ScalingMethod | str = 'linear',
                      show: bool = True,
                      new_figure: bool = True,
                      use_rust: bool = True,
                      data: ArrayLike | None = None,
                      ignore_imaginary: bool = False):
        """Show the powder spectrum.

        Parameters
        ----------
        path : Path1D or EmpiricalPath1D or array_like or tuple
            1D path or array of ``|q|`` magnitudes in inverse Angstroms.
        n_samples : int, default=100
            Approximate number of spherical orientation samples per ``|q|`` shell.
        method : SphericalPointGeneratorType or str, default='fibonacci'
            Spherical point sampling method.
        min_energy : float, optional
            Minimum energy for histogramming in meV.
        max_energy : float, optional
            Maximum energy for histogramming in meV.
        n_energy_bins : int, optional
            Number of energy bins.
        energy_stddev : float, optional
            Standard deviation for Gaussian energy broadening in meV.
        random_seed : int, optional
            Random seed for stochastic spherical sampling methods.
        scaling_method : ScalingMethod or str, default='linear'
            Intensity scaling method (``'linear'`` or ``'log'``).
        show : bool, default=True
            If ``True``, call ``matplotlib.pyplot.show()``.
        new_figure : bool, default=True
            If ``True``, create a new matplotlib figure.
        use_rust : bool, default=True
            Whether to use the Rust backend for calculation.
        data : array_like, optional
            Precomputed 2D spectrum array to plot instead of recomputing.
        ignore_imaginary : bool, default=False
            If ``True``, zero out intensities of modes with significant imaginary
            energy components.

        Returns
        -------
        matplotlib.figure.Figure or None
            The matplotlib figure if ``show=False``, otherwise ``None``.
        """
        if not isinstance(path, Path1DBase):
            path = EmpiricalPath1D(path)

        if data is None or min_energy is None or max_energy is None:
            q, e, data = self.spectrum(path, n_samples, method,
                                       min_energy, max_energy, n_energy_bins, energy_stddev,
                                       random_seed, use_rust, ignore_imaginary)
        else:
            q, e = path.q_values(), np.linspace(min_energy, max_energy, data.shape[1])

        if isinstance(scaling_method, str):
            scaling_method = ScalingMethod(scaling_method)

        match scaling_method:
            case ScalingMethod.LOG:
                zeros = data <= 0
                min_value = np.log10(np.min(data[~zeros][:]))

                data = np.log10(data)
                data[zeros] =  min_value - 1


        import matplotlib.pyplot as plt

        if new_figure:
            if plt.fignum_exists("Powder Spectrum"):
                fig, ax = plt.subplots()
            else:
                fig, ax = plt.subplots(num="Powder Spectrum")
        else:
            fig, ax = plt.gcf(), plt.gca()

        ax.pcolormesh(q, e, data.T, shading='nearest')
        ax.set_xlabel(r"$|q| (\AA^{-1})$")
        ax.set_ylabel("Energy (meV)")
        ax.set_aspect(0.2)

        if show:
            plt.show()
        else:
            return fig

    def _sample_serialise(self, context: SPWSerialisationContext):
        return {
            "hamiltonian": self.hamiltonian._serialise(context)
        }

    @staticmethod
    @expects_keys("hamiltonian")
    def _sample_deserialise(json, context: SPWDeserialisationContext):
        return Powder(Hamiltonian._deserialise(json["hamiltonian"], context))

class MagneticFieldPowder(Sample1D):
    """Powder in magnetic field (slightly different parallelisation method to Powder)."""

    def __init__(self, hamiltonian: Hamiltonian):
        raise NotImplementedError("Powder in magnetic field not implemented yet")

class TwoMagnon(Sample3D):
    """Two magnon excitations in a single crystal."""

    def __init__(self, hamiltonian: Hamiltonian):
        raise NotImplementedError("Two magnon calculation not implemented yet")


_all_samples = [SingleCrystal, Twin, Multidomain, Powder]
_all_samples_lookup = {cls.sample_name: cls for cls in _all_samples}
