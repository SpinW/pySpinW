"""Different Kinds of samples"""

from abc import ABC, abstractmethod
from enum import Enum
from typing import Sequence, Callable

import numpy as np
from numpy._typing import ArrayLike

import matplotlib.pyplot as plt

from pyspinw.units import IntensityUnits
from pyspinw.calculations.spherical_integration import SphericalPointGeneratorType, point_generator
from pyspinw.checks import check_sizes
from pyspinw.hamiltonian import Hamiltonian, ParametrizationType, omegasum, egrid
from pyspinw.path import Path, Path1D, Path1DBase, EmpiricalPath1D
from pyspinw.tolerances import tolerances


# pylint: disable=R0903

class Sample(ABC):
    """Representation of the macrostructure of a sample used in an experiment (Twin, Powder etc)"""

    def __init__(self, hamiltonian: Hamiltonian):
        self.hamiltonian = hamiltonian

    def _energies(self,
                 path: Path,
                 field: ArrayLike | None = None,
                 use_rust: bool=True):

        return self.hamiltonian.sorted_positive_energies(path.q_points(), field=field, use_rust=use_rust)

    @abstractmethod
    def _energies_and_intensities(self,
                                 points: np.ndarray,
                                 field: ArrayLike | None = None,
                                 use_rust: bool = True,
                                 use_rotating: bool=True,
                                 intensity_unit: IntensityUnits | str='cell',):
        """ Abstract method to get the energies and intensities """

class Sample3D(Sample):
    """ Sample where the direction of q matters"""

    @check_sizes(field=(3,), allow_nones=True)
    def energies(self,
                 path: Path,
                 field: ArrayLike | None = None,
                 use_rust: bool = True) -> list[np.ndarray]:
        """ Get the energies along a specified path """
        return self._energies(path, field, use_rust)


    def energies_and_intensities(self,
                                 path: Path,
                                 field: ArrayLike | None = None,
                                 use_rust: bool = True,
                                 use_rotating: bool = True,
                                 intensity_unit: IntensityUnits | str = 'cell'):

        """ Get energy and intensity data"""

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
                       intensity_unit: IntensityUnits | str = 'cell',
                       scale: str='linear'):
        """ Create a spaghetti diagram with intensity as colorfill overplotted by mode energies """
        if new_figure:
            fig, ax = plt.subplots()
        else:
            fig = plt.gcf()
            ax = fig.get_axes()[0]

        x_values = path.x_values()
        energy, intensity = omegasum(*self._energies_and_intensities(path.q_points(), field,
                                                                    use_rust, use_rotating, intensity_unit))
        if evect is None:
            emax = np.nanmax(np.real(energy))
            denom = 10**np.floor(np.log10(emax))
            evect = np.linspace(0, np.ceil(emax / denom) * denom, 100)
        if dE is None:
            dE = np.mean(np.diff(evect)) * 2.35
        spec = egrid(energy, intensity, evect, dE)
        if vmax is None:
            vmax = np.nanmax(spec[:, np.where(evect > max(np.max(evect)/10, 0.1))]) / 10.
        mesh = ax.pcolormesh(x_values, evect, spec.T, vmin=vmin, vmax=vmax)
        ax.plot(x_values, np.real(energy), 'k')
        if (np.imag(energy) > 0.01).any():
            ax.plot(x_values, np.imag(energy), 'or')
        ax.set_ylim(0, np.max(evect))
        fig.colorbar(mesh, ax=ax)
        path.format_plot(ax)
        ax.set_ylabel('Magnon Energy (meV)')

        if show:
            plt.show()
        else:
            return fig


class Sample1D(Sample):
    """ Sample where the direction of q does not matter

    An important difference from Sample3D is there's no option to just get energies,
    this is because they will form a continuum in 1D. The 3D q-energy curves can
    be obtained from the Hamiltonian, they're just not exposed at this level
    """



class CrystalDomain:
    """ Orientation and amount of crystalline domain (e.g. one of a twin, or a domain in the traditional sense) """

    @check_sizes(transformation=(3,3), force_numpy=True)
    def __init__(self, transformation: ArrayLike, weighting: float):

        # check the transformation is valid
        if not np.allclose(transformation @ transformation.T, np.eye(3), atol=tolerances.IS_ZERO_TOL):
            raise ValueError(f"Expected crystal domain transformation to be orthonormal, got {transformation}")

        if weighting < 0:
            raise ValueError(f"Weightings must be positive, got {weighting}")

        self.transformation = transformation
        self.weighting = weighting

CrystalDomainLike = CrystalDomain | tuple[ArrayLike, float]


def _domain_like_to_domain(domain_like: CrystalDomainLike) -> CrystalDomain:
    """ Converts CrystalDomainLike to CrystalDomain """
    if isinstance(domain_like, CrystalDomain):
        return domain_like
    else:
        return CrystalDomain(*domain_like)


class SingleCrystal(Sample3D):
    """Specifies a single crystal sample"""

    def __init__(self, hamiltonian: Hamiltonian):

        super().__init__(hamiltonian=hamiltonian)

    def _energies_and_intensities(self,
                                 q_points: np.ndarray,
                                 field: ArrayLike | None = None,
                                 use_rust: bool = True,
                                 use_rotating: bool = True,
                                 intensity_unit: IntensityUnits | str = 'cell'):

        """ Get energy and intensity data"""

        # Note: we don't use the _ method because we want to ignore non-magnetic sites
        return self.hamiltonian.energies_and_intensities(
            q_points, field=field, use_rust=use_rust, use_rotating=use_rotating, intensity_unit=intensity_unit)


class Multidomain(Sample3D):
    """Sample consisting of multiple domains"""

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

            transformed_q = transformation @ q_vectors
            transformed_field = None if field is None else transformation @ field

            output += self.hamiltonian.sorted_positive_energies(
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

            energy, intensities = self.hamiltonian.energies_and_intensities(
                        transformed_q,
                        field=transformed_field,
                        use_rust=use_rust)

            output_energies.append(energy)
            for i in range(n_q):
                output_intensities[i] += [weight * intensity for intensity in intensities[i]]

        output_energies = np.concatenate(output_energies, axis=1)

        return output_energies, output_intensities



class Twin(Multidomain):
    """Specify a twinned crystal.

    Special case of multidomain
    """

    def __init__(self,
                 hamiltonian: Hamiltonian,
                 relative_rotation: ArrayLike,
                 second_twin_fraction: float):

        if not (0 <= second_twin_fraction <= 1):
            raise ValueError("Expected weighting to be between 0 and 1")

        first_twin_fraction = 1 - second_twin_fraction

        super().__init__(
            hamiltonian=hamiltonian,
            domains=[CrystalDomain(np.eye(3), first_twin_fraction),
                     CrystalDomain(relative_rotation, second_twin_fraction)])

class ScalingMethod(Enum):
    """ Scaling methods for plots"""

    LINEAR = 'linear'
    LOG = 'log'

class ParameterizedPowderSpectrum:
    """ A powder spectrum that has been parameterised"""

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
                 find_ground_state_with: dict | None=None):

        self._powder = Powder
        self._parameters = parameters

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

        self.parameterised_hamiltonian = powder.hamiltonian.parameterize(
                                *parameters,
                                           find_ground_state_with=find_ground_state_with)

    def __call__(self, *parameters):
        """ Evaluate the powder spectrum at specified parameter values """
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
            use_rust=self._use_rust)[2]


class Powder(Sample1D):
    """Sample is a powder"""

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
                 use_rust: bool = True):
        """ Get the powder spectrum """
        if not isinstance(path, Path1DBase):
            path = EmpiricalPath1D(path)

        generator = point_generator(method)(n_samples, seed = random_seed)

        chunk_size = generator.actual_n_points
        qs = path.q_values()
        points = np.concatenate([generator.points*q for q in qs], axis=0)

        energies, intensities = self.hamiltonian.energies_and_intensities(points, use_rust=use_rust)

        energies = np.real(np.array(energies))
        intensities = np.real(np.array(intensities))

        positive_energies = energies > 0


        # Historgramming parameters
        min_energy = np.min(energies[positive_energies]) if min_energy is None else min_energy
        max_energy = np.max(energies[positive_energies]) if max_energy is None else max_energy
        n_energy_bins = n_samples // 4 if n_energy_bins is None else n_energy_bins

        energy_bin_edges = np.linspace(min_energy, max_energy, n_energy_bins + 1)
        energy_bin_centres = 0.5 * (energy_bin_edges[1:] + energy_bin_edges[:-1])

        # Output map
        output = np.zeros((path.n_points, n_energy_bins))

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


            output[i, :] = binned

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
            find_ground_state_with: dict | None = None):
        """ Get the powder spectrum as a function of parameters """
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
                        find_ground_state_with=find_ground_state_with)

    def show_spectrum(self,
                      path: Path1D | EmpiricalPath1D | ArrayLike,
                      n_samples: int = 100,
                      method: SphericalPointGeneratorType | str = "fibonacci",
                      min_energy: float | None = None,
                      max_energy: float | None = None,
                      n_energy_bins: int | None = None,
                      energy_stddev: float | None = None,
                      random_seed: int | None = None,
                      scaling_method: ScalingMethod | str = 'linear',
                      show_plot: bool = True,
                      new_figure: bool = True,
                      use_rust: bool = True):
        """ Show the powder spectrum """
        if not isinstance(path, Path1DBase):
            path = EmpiricalPath1D(path)

        q, e, data = self.spectrum(path, n_samples, method,
                                   min_energy, max_energy, n_energy_bins, energy_stddev,
                                   random_seed, use_rust)

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
            plt.figure("Powder Spectrum")

        ax = plt.gca()

        plt.imshow(data.T[::-1, :], extent=(q[0], q[-1], e[0], e[-1]))
        ax.set_aspect(0.2)

        if show_plot:
            plt.show()



class MagneticFieldPowder(Sample1D):
    """ Powder in magnetic field (slightly different parallelisation method to Powder)"""

class TwoMagnon(Sample3D):
    """ Two magnon excitations in a single crystal"""
