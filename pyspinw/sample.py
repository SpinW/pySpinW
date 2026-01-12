"""Different Kinds of samples"""

from abc import ABC, abstractmethod

import numpy as np
from numpy._typing import ArrayLike

from pyspinw.calculations.spherical_integration import SphericalPointGeneratorType, point_generator
from pyspinw.checks import check_sizes
from pyspinw.hamiltonian import Hamiltonian
from pyspinw.path import Path, Path1D
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


class Sample3D(Sample):
    """ Sample where the direction of q matters"""

    @check_sizes(field=(3,), allow_nones=True)
    def energies(self,
                 path: Path,
                 field: ArrayLike | None = None,
                 use_rust: bool = True) -> list[np.ndarray]:
        """ Get the energies along a specified path """

        return self._energies(path, field, use_rust)

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



class Powder(Sample1D):
    """Sample is a powder"""

    def __init__(self, hamiltonian: Hamiltonian):

        super().__init__(hamiltonian)

    def spectrum(self,
                 path: Path1D,
                 n_samples: int=100,
                 method: SphericalPointGeneratorType | str = "fibonacci",
                 min_energy: float | None = None,
                 max_energy: float | None = None,
                 n_energy_bins: int | None = None,
                 random_seed: int | None = None,
                 use_rust: bool = True):

        """ Get the powder spectrum """

        generator = point_generator(method)(n_samples, seed = random_seed)

        chunk_size = generator.actual_n_points
        qs = path.q_values()
        points = np.concatenate([generator.points*q for q in qs], axis=0)

        energies, intensities = self.hamiltonian.energies_and_intensities(points)

        positive_energies = energies > 0


        # Historgramming parameters
        min_energy = np.min(energies[positive_energies]) if min_energy is None else min_energy
        max_energy = np.max(energies[positive_energies]) if max_energy is None else max_energy
        n_energy_bins = n_samples // 4 if n_energy_bins is None else n_energy_bins

        energy_bin_edges = np.linspace(min_energy, max_energy, n_energy_bins + 1)
        energy_bin_centres = 0.5 * (energy_bin_edges[1:] + energy_bin_edges[:-1])

        # Output map
        output = np.zeros((path.resolution, n_energy_bins))

        # at each q, bin the energies with weights of the intensities, but ignore the negative ones
        for i, q in enumerate(path.q_values()):

            # Get the part of the data that's relevant
            energy = energies[i*chunk_size:(i+1)*chunk_size]
            intensity = intensities[i*chunk_size:(i+1)*chunk_size]
            positives = positive_energies[i*chunk_size:(i+1)*chunk_size]

            energy = energy[positives]
            intensity = intensity[positives]

            # Do the binning
            binned, bin_edges = np.histogram(energy, bins=energy_bin_edges, weights=intensity)

            output[i, :] = binned

        return path.q_values(), energy_bin_centres, output

    def show_spectrum(self,
                      path: Path1D,
                      n_samples: int = 100,
                      method: SphericalPointGeneratorType | str = "fibonacci",
                      min_energy: float | None = None,
                      max_energy: float | None = None,
                      n_energy_bins: int | None = None,
                      random_seed: int | None = None,
                      show_plot: bool = True,
                      new_figure: bool = True,
                      use_rust: bool = True):

        q, e, data = self.spectrum(path, n_samples, method, min_energy, max_energy, n_energy_bins, random_seed, use_rust)

        import matplotlib.pyplot as plt

        if new_figure:
            plt.figure("Powder Spectrum")

        ax = plt.gca()

        plt.imshow(data, extent=(q[0], q[-1], e[0], e[-1]))
        ax.set_aspect(0.2)

        if show_plot:
            plt.show()



class MagneticFieldPowder(Sample1D):
    pass

class TwoMagnon(Sample3D):
    """ Two magnon excitations in a single crystal"""
